"""RM16's second question: is a frame kept beyond what its host gene's conservation gives for free?

Three nulls, each per peptide, reading and class (a class with at least SPECIES_MIN assessed species):

window      H4's null made exact. The frames of the same length in the same gene, on the same strand and, inside a
            protein gene, in the same frame, which are ORFs in human under the same reading (statistic "kept": a
            start codon and no stop in n codons) or merely stop-free in human (statistic "open", H4's own), every
            one of them except those overlapping the peptide in its frame. The statistic is the share of the class'
            assessed species where the window is kept (open). p = (1 + windows at or above the observed share) /
            (1 + windows). Fewer than WINDOWS_MIN windows: not assessed.
synonymous  protein genes only. Each species' host codons are shuffled among synonymous positions (shuffle.py),
            REPLICATES times; the frame is read in each shuffle at the species' own window. The statistic is the
            share of species kept; p = (1 + replicates at or above the observed share) / (1 + replicates).
dinucleotide rRNA genes only. Each species' window (its first 3n bases from the start column) is shuffled keeping
            its dinucleotide counts, REPLICATES times; the same statistic and p.
The shuffle nulls draw each species independently, while related species inherit their frames: at the level of
species they overstate the evidence. They are therefore also computed with one species per family (the family's
assessed species with the smallest SHA-256 of its slug), and that level is the one read.
q: Benjamini-Hochberg over peptides and classes, within null, statistic, level and reading.
Writes nulls.tsv.
"""
import hashlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
from statsmodels.stats.multitest import multipletests

from . import codes, frames, shuffle
from .catalogue import PEPTIDES, Peptide

REPLICATES = 1000
WINDOWS_MIN = 20
SEED = 20261115


def window_starts(pc: frames.Placed, reading: str, human_code: int, statistic: str) -> list[int]:
    """Start positions (strand-oriented human sequence) of the null windows."""
    host, h = pc.host, pc.host.row("homo-sapiens")
    L, n = int(host.length[h]), pc.n
    seq = host.sequence(h)
    table = codes.tables(human_code, codes.species_reading(reading, pc.pep.reading))
    protein = pc.pep.alignment.endswith(".codon")
    out = []
    for q in range(0, L - 3 * n + 1):
        if protein and q % 3 != pc.q0 % 3:
            continue
        if (q - pc.q0) % 3 == 0 and abs(q - pc.q0) < 3 * n:      # overlaps the peptide in its own frame
            continue
        aa = codes.translate(seq[q:q + 3 * n], table)
        if "*" in aa or "X" in aa:
            continue
        if statistic == "kept" and not table.is_start(seq[q:q + 3]):
            continue
        out.append(q)
    return out


def class_shares(state: np.ndarray, open_: np.ndarray, cls: np.ndarray, classes: list[str], statistic: str) -> dict:
    """{class: (assessed, share)} for one evaluation of the frame in every species."""
    hit = (state == frames.KEPT) if statistic == "kept" else open_
    assessed = state != frames.NOT_ASSESSED
    out = {}
    for c in classes:
        sel = assessed & (cls == c)
        k = int(sel.sum())
        out[c] = (k, float(hit[sel].mean()) if k else float("nan"))
    return out


def window_null(pc: frames.Placed, reading: str, code_of: dict[str, int], cls: np.ndarray, classes: list[str]) -> list[dict]:
    tabs = frames.tables_for(pc, reading, code_of)
    rd = frames.reader(pc.host, tabs)
    obs = frames.state_at(pc.host, rd, pc.c0, pc.cols(), pc.n)
    rows = []
    for statistic in ("kept", "open"):
        starts = window_starts(pc, reading, code_of["homo-sapiens"], statistic)
        o = class_shares(obs["state"], obs["open"], cls, classes, statistic)
        null = {c: [] for c in classes}
        for q in starts:
            st = frames.state_at(pc.host, rd, int(pc.hcols[q]), pc.cols(q), pc.n)
            for c, (k, s) in class_shares(st["state"], st["open"], cls, classes, statistic).items():
                if k >= frames.SPECIES_MIN:
                    null[c].append(s)
        for c in classes:
            k, s = o[c]
            v = np.array(null[c])
            base = {"peptide": pc.pep.name, "reading": reading, "null": "window", "statistic": statistic,
                    "level": "species", "class": c, "n": k, "observed": s, "null_size": len(v)}
            if k < frames.SPECIES_MIN:
                rows.append({**base, "null_mean": None, "p": None, "status": f"not assessed: under {frames.SPECIES_MIN} assessed species"})
            elif len(v) < WINDOWS_MIN:
                rows.append({**base, "null_mean": float(v.mean()) if len(v) else None, "p": None,
                             "status": f"not assessed: {len(v)} null windows, under {WINDOWS_MIN}"})
            else:
                rows.append({**base, "null_mean": float(v.mean()), "p": float((1 + (v >= s - 1e-12).sum()) / (1 + len(v))),
                             "status": "tested"})
    return rows


def representatives(slugs: list[str], meta: dict) -> set[str]:
    """One species per family: the one with the smallest SHA-256 of its slug."""
    best: dict[str, tuple[str, str]] = {}
    for s in slugs:
        unit = meta[s]["family"] or meta[s]["genus"] or s
        h = hashlib.sha256(s.encode()).hexdigest()
        if unit not in best or h < best[unit][0]:
            best[unit] = (h, s)
    return {s for _, s in best.values()}


def shuffle_draws(pc: frames.Placed, code_of: dict[str, int], replicates: int, seed: int) -> dict:
    """Per species (assessed under the own reading's observation of the window position), the kept and open flags of
    the unshuffled frame and of `replicates` shuffles, for both readings. -> {"slugs": [...], reading: {"kept_obs",
    "open_obs", "kept", "open"}} with arrays (species,) and (replicates, species); species the shuffle cannot read
    (see `why`) are listed in "skipped" with the reason."""
    host, pep, n = pc.host, pc.pep, pc.n
    m = frames.needed(n)
    own_tabs = frames.tables_for(pc, "own", code_of)
    obs = frames.state_at(host, frames.reader(host, own_tabs), pc.c0, pc.cols(), n)
    protein = pep.alignment.endswith(".codon")
    slugs, skipped = [], []
    res = {r: {"kept_obs": [], "open_obs": [], "kept": [], "open": []} for r in frames.READINGS}
    for i, slug in enumerate(host.slugs):
        if obs["state"][i] == frames.NOT_ASSESSED:
            continue
        p, L = int(obs["p"][i]), int(host.length[i])
        rng = shuffle.rng_for(seed, pep.name, slug)
        bases = host.nt[i, :L]
        if protein:
            if L % 3:
                skipped.append((slug, "host length not a multiple of 3")); continue
            cod = codes.codon_index(bases)[0::3][:L // 3]
            k0, k1 = p // 3, (p + 3 * m - 1) // 3
            if k1 >= len(cod) or (cod[k0:k1 + 1] < 0).any():
                skipped.append((slug, "unknown base in the window's host codons")); continue
            host_table = codes.tables(code_of[slug], "own")       # synonymy under the species' own code
            aa = shuffle.residues(cod, host_table)
            draws = shuffle.synonymous(cod, aa, replicates, rng)[:, k0:k1 + 1]
            real = cod[k0:k1 + 1][None, :]
            off = p - 3 * k0
            alt_draws = shuffle.frame_codons(shuffle.codons_to_bases(draws), off, m)
            alt_real = shuffle.frame_codons(shuffle.codons_to_bases(real), off, m)
        else:
            seg = bases[p:p + 3 * min(n, (L - p) // 3)]
            if len(seg) < 3 * m or (seg > 3).any():
                skipped.append((slug, "window shorter than the needed codons or with an unknown base")); continue
            s = "".join("ACGT"[b] for b in seg)
            alt_draws = shuffle.frame_codons(shuffle.dinucleotide_many(s, replicates, rng), 0, m)
            alt_real = shuffle.frame_codons(seg[None, :], 0, m)
        slugs.append(slug)
        for r in frames.READINGS:
            t = codes.tables(code_of[slug], codes.species_reading(r, pep.reading))
            k, o = shuffle.kept_open(alt_draws, t)
            kr, orr = shuffle.kept_open(alt_real, t)
            res[r]["kept"].append(k); res[r]["open"].append(o)
            res[r]["kept_obs"].append(bool(kr[0])); res[r]["open_obs"].append(bool(orr[0]))
    out = {"slugs": slugs, "skipped": skipped}
    for r in frames.READINGS:
        out[r] = {"kept_obs": np.array(res[r]["kept_obs"], dtype=bool), "open_obs": np.array(res[r]["open_obs"], dtype=bool),
                  "kept": np.stack(res[r]["kept"], axis=1) if slugs else np.zeros((replicates, 0), bool),
                  "open": np.stack(res[r]["open"], axis=1) if slugs else np.zeros((replicates, 0), bool)}
    return out


def shuffle_rows(pep: Peptide, draws: dict, meta: dict, classes: list[str], replicates: int) -> list[dict]:
    kind = "synonymous" if pep.alignment.endswith(".codon") else "dinucleotide"
    slugs = np.array(draws["slugs"], dtype=object)
    cls = np.array([meta[s]["class"] for s in slugs], dtype=object)
    rows = []
    for reading in frames.READINGS:
        d = draws[reading]
        for statistic in ("kept", "open"):
            obs, null = d[f"{statistic}_obs"], d[statistic]
            for level in ("species", "family"):
                for c in classes:
                    sel = cls == c
                    if level == "family":
                        reps = representatives(list(slugs[sel]), meta)
                        sel = sel & np.isin(slugs, list(reps))
                    k = int(sel.sum())
                    base = {"peptide": pep.name, "reading": reading, "null": kind, "statistic": statistic,
                            "level": level, "class": c, "n": k, "null_size": replicates}
                    if k < frames.SPECIES_MIN:
                        rows.append({**base, "observed": float(obs[sel].mean()) if k else None, "null_mean": None,
                                     "p": None, "status": f"not assessed: under {frames.SPECIES_MIN} species"})
                        continue
                    o = float(obs[sel].mean())
                    v = null[:, sel].mean(1)
                    rows.append({**base, "observed": o, "null_mean": float(v.mean()),
                                 "p": float((1 + (v >= o - 1e-12).sum()) / (1 + replicates)), "status": "tested"})
    return rows


def _task(args):
    data, name, replicates, seed = args
    sp = frames.species_table(data)
    meta = {r["slug"]: r for r in sp.iter_rows(named=True)}
    code_of = {s: int(r["transl_table"]) for s, r in meta.items()}
    pep = next(p for p in PEPTIDES if p.name == name)
    pc = frames.placed(data, pep, set(meta))
    cls = np.array([meta[s]["class"] for s in pc.host.slugs], dtype=object)
    classes = sorted(set(cls))
    rows = []
    for reading in frames.READINGS:
        rows += window_null(pc, reading, code_of, cls, classes)
    draws = shuffle_draws(pc, code_of, replicates, seed)
    rows += shuffle_rows(pep, draws, meta, classes, replicates)
    reasons: dict[str, int] = {}
    for _, why in draws["skipped"]:
        reasons[why] = reasons.get(why, 0) + 1
    note = "; ".join(f"{why}: {k}" for why, k in sorted(reasons.items()))
    return name, rows, note


def main(data: Path, out: Path, replicates: int = REPLICATES, seed: int = SEED, jobs: int = 1,
         peptides: list[str] | None = None) -> pl.DataFrame:
    out.mkdir(parents=True, exist_ok=True)
    names = [p.name for p in PEPTIDES if peptides is None or p.name in peptides]
    tasks = [(data, nm, replicates, seed) for nm in names]
    if jobs > 1:
        with ProcessPoolExecutor(jobs) as ex:
            results = list(ex.map(_task, tasks))
    else:
        results = [_task(t) for t in tasks]
    rows = []
    for name, r, note in results:
        rows += r
        print(f"{name}: shuffle null skipped species: {note or 'none'}")
    res = pl.DataFrame(rows, infer_schema_length=None).with_columns(pl.lit(seed).alias("seed"))
    res = bh(res)
    res = res.sort(["null", "statistic", "level", "reading", "peptide", "class"])
    res.write_csv(out / "nulls.tsv", separator="\t", float_precision=6)
    return res


def bh(res: pl.DataFrame) -> pl.DataFrame:
    """q within null, statistic, level and reading, over the tested rows."""
    q = np.full(res.height, np.nan)
    keys = ["null", "statistic", "level", "reading"]
    pdf = res.with_row_index("_i")
    for _, g in pdf.filter(pl.col("status") == "tested").group_by(keys, maintain_order=True):
        q[g["_i"].to_numpy()] = multipletests(g["p"].to_numpy(), method="fdr_bh")[1]
    return res.with_columns(pl.Series("q", q).fill_nan(None))
