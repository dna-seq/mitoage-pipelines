"""H31: the cybrids and the human variants scored from a masked DNA model's likelihood table of the host's genome.

The table (`slurm/h31-kron-table.py`, KronGenome) holds, for each place of a host's mitochondrial genome, the
log-probability of each base there when it is hidden. A change of the host's base to another is scored as the
log-probability of the host's base minus that of the other: positive when the other is the less expected.

The rules are in docs/HYPOTHESES.md (H31) and were written before this code.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner

from . import cybrid, cybrid_dna, variant_stab

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
PRIMARY = ("kron_cds", "kron_cds_mean")
PAIR = ("kron_rrna", "kron_trna", "kron_dloop")
BASELINE = ("d_all", "f_all", "nt_cds_d", "nt_cds_f")
T3_MIN = 5
AUC_MIN = 0.70
KINDS = ("protein", "tRNA", "rRNA", "control region", "other")


def read_table(path: Path) -> tuple[str, np.ndarray]:
    """(the genome as the table has it, [n, 4] log-probabilities of A, C, G, T; NaN where a place was not scored)."""
    t = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""])
    if list(t.pos) != list(range(1, len(t) + 1)):
        raise ValueError(f"{path}: the places are not 1..n in order")
    return "".join(t.ref), t[[f"lp_{b}" for b in BASES]].to_numpy(float)


def drop(lp: np.ndarray, place: int, ref: str, alt: str) -> float:
    """The log-probability of `ref` minus that of `alt` at the 0-based place; NaN when either is not A, C, G or T or
    the place was not scored."""
    if ref not in BASES or alt not in BASES:
        return float("nan")
    return float(lp[place, BASES.index(ref)] - lp[place, BASES.index(alt)])


def stretches(gb: Path, accession: str) -> dict[tuple, tuple[str, np.ndarray, int]]:
    """{("cds", gene) | ("rrna", "12S"|"16S") | ("trna", (product, i)) | ("dloop",): (the stretch read 5' to 3' on
    its own strand, the record's 0-based place of each of its bases, its strand)} of one record. The control region
    is read on the plus strand in the genome's order."""
    rec = next((r for r in SeqIO.parse(str(gb), "genbank") if r.id == accession), None)
    if rec is None:
        raise SystemExit(f"{accession} is not in {gb}: `mitoage download cybrid`")
    out, seen = {}, {}
    for f in rec.features:
        product = f.qualifiers.get("product", [None])[0]
        if f.type == "CDS":
            g = cybrid.gene_name(f.qualifiers.get("gene", [None])[0], product)
            key = ("cds", g) if g else None
        elif f.type == "rRNA":
            key = ("rrna", cybrid_dna.RRNA.get((product or "").strip().lower()))
        elif f.type == "tRNA":
            i = seen.get(product, 0); seen[product] = i + 1
            key = ("trna", (product, i))
        elif f.type == "D-loop" or (f.type == "misc_feature"
                                    and f.qualifiers.get("note", [""])[0].strip().lower() == "control region"):
            parts = sorted(f.location.parts, key=lambda q: int(q.start) == 0)
            places = np.concatenate([np.arange(int(q.start), int(q.end)) for q in parts])
            out[("dloop",)] = ("".join(str(rec.seq[i]) for i in places).upper(), places, 1)
            continue
        else:
            continue
        if key is None:
            continue
        loc = cybrid_dna.trna_location(rec.id, f) if f.type == "tRNA" else f.location
        places = np.array(list(loc), dtype=np.int64)                 # in the stretch's own 5' to 3' order
        out[key] = (str(loc.extract(rec.seq)).upper(), places, int(loc.strand))
    if ("dloop",) not in out and "dloop" in cybrid_dna.RECORD_FIXES.get(rec.id, {}):
        s, e = cybrid_dna.RECORD_FIXES[rec.id]["dloop"]
        out[("dloop",)] = (str(rec.seq[s:e]).upper(), np.arange(s, e), 1)
    return out


def stretch_score(host: tuple[str, np.ndarray, int], donor: str, genome: str, lp: np.ndarray,
                  al: PairwiseAligner) -> tuple[float, int, int]:
    """(the sum of the drops over the columns where host and donor both hold a base and differ, how many such columns
    were scored, the columns with a gap) for one stretch. The table is of the plus strand: on the minus strand both
    bases are complemented. Refuses a host base that is not the table's."""
    seq, places, strand = host
    total, n, gaps, i = 0.0, 0, 0, -1
    for x, y in cybrid_dna.columns(seq, donor, al):
        if x != "-":
            i += 1
        if x == "-" or y == "-":
            gaps += 1
            continue
        if x == y:
            continue
        ref, alt = (x, y) if strand == 1 else (x.translate(COMPLEMENT), y.translate(COMPLEMENT))
        if genome[places[i]] != ref:
            raise ValueError(f"place {places[i] + 1} is {genome[places[i]]} in the table, not {ref}")
        d = drop(lp, int(places[i]), ref, alt)
        if not np.isnan(d):
            total += d; n += 1
    return total, n, gaps


def pair_rows(gb: Path, host: str, donor: str, genome: str, lp: np.ndarray, al: PairwiseAligner | None = None) -> list[dict]:
    al = al or cybrid_dna.aligner()
    h = stretches(gb, cybrid.GENOMES[host])
    d = stretches(gb, cybrid.GENOMES[donor])
    own = {g: stretches(gb, a)[("cds", g)] for (s, g), a in cybrid.OWN_LINE.items() if s == donor}

    def over(keys):
        got = [stretch_score(h[k], (own[k[1]] if k[0] == "cds" and k[1] in own else d[k])[0], genome, lp, al) for k in keys]
        return sum(x for x, _, _ in got), sum(n for _, n, _ in got), sum(g for _, _, g in got)

    pair = {}
    for name, keys in (("kron_rrna", [("rrna", "12S"), ("rrna", "16S")]), ("kron_trna", sorted(k for k in h if k[0] == "trna"))):
        if not set(keys) <= set(d):
            raise ValueError(f"{host} / {donor}: the genes of {name} do not pair")
        crossed = sorted(str(k[1]) for k in keys if h[k][2] != d[k][2])
        if crossed:
            raise ValueError(f"{host} / {donor}: {crossed} are on different strands in the two records")
        s, n, _ = over(keys)
        pair[name], pair[name + "_mean"] = s, (s / n if n else np.nan)
    if ("dloop",) in h and ("dloop",) in d:
        s, n, g = over([("dloop",)])
        pair.update(kron_dloop=s, kron_dloop_mean=s / n if n else np.nan, kron_dloop_gaps=g)
    else:
        pair.update(kron_dloop=np.nan, kron_dloop_mean=np.nan, kron_dloop_gaps=np.nan)
    rows = []
    for cx, genes in cybrid.COMPLEX_GENES.items():
        s, n, g = over([("cds", x) for x in genes])
        rows.append({"host": host, "donor": donor, "complex": cx, "kron_cds": s, "kron_cds_mean": s / n if n else np.nan,
                     "kron_cds_n": n, "kron_cds_gaps": g, **pair})
    return rows


def verdict(summary: pd.DataFrame) -> pd.DataFrame:
    """The cybrid rule: each primary score against the best of the baselines on T1 (lower), T2 (higher) and T3 (more
    pairs named, and at least T3_MIN)."""
    def value(test, s):
        v = summary[(summary.test == test) & (summary.score == s)].value
        return float(v.sum()) if test == "T3" else float(v.iloc[0])
    bar = {"T1": min(value("T1", b) for b in BASELINE), "T2": max(value("T2", b) for b in BASELINE),
           "T3": max(T3_MIN - 1, max(value("T3", b) for b in BASELINE))}
    rows = []
    for s in PRIMARY:
        t1, t2, t3 = value("T1", s), value("T2", s), value("T3", s)
        beats = [t1 < bar["T1"], t2 > bar["T2"], t3 > bar["T3"]]
        rows.append({"score": s, "T1": t1, "T1_bar": bar["T1"], "T2": t2, "T2_bar": bar["T2"], "T3": t3, "T3_bar": bar["T3"],
                     "verdict": "counts" if all(beats) else "partial, claimed as nothing" if sum(beats) == 2 else "does not count"})
    return pd.DataFrame(rows)


def cybrids(gb: Path, tables: Path, out: Path, arm: str = "mito", measurements: Path = cybrid.MEASUREMENTS) -> pd.DataFrame:
    """kron_scores.tsv, h31_calibration.tsv, h31_summary.tsv and h31_verdict.tsv in `out`, which must hold scores.tsv
    (`mitoport cybrid`) and dna_scores.tsv (`mitoport cybrid-dna`) for the baselines."""
    key = ["host", "donor", "complex"]
    for f in ("scores.tsv", "dna_scores.tsv"):
        if not (out / f).is_file():
            raise SystemExit(f"{out / f} is missing: run `mitoage mitoport cybrid` and `mitoage mitoport cybrid-dna` first")
    meas = cybrid.read_measurements(measurements)
    al, rows = cybrid_dna.aligner(), []
    for host in sorted(set(meas.host)):
        genome, lp = read_table(tables / f"{cybrid.GENOMES[host]}.{arm}.tsv")
        for donor in sorted(set(meas.donor[meas.host == host])):
            rows += pair_rows(gb, host, donor, genome, lp, al)
    kron = pd.DataFrame(rows)
    kron.to_csv(out / "kron_scores.tsv", sep="\t", index=False, float_format="%.4f")
    scores = (pd.read_csv(out / "scores.tsv", sep="\t")[key + ["d_all", "f_all"]]
              .merge(pd.read_csv(out / "dna_scores.tsv", sep="\t")[key + ["nt_cds_d", "nt_cds_f"]], on=key)
              .merge(kron, on=key, how="left"))
    names = PRIMARY + PAIR + BASELINE
    cal, summary = cybrid.summarise(scores, meas, names, t2b=PRIMARY + BASELINE, t3=PRIMARY + BASELINE,
                                    t4=("kron_cds", *PAIR, "d_all", "nt_cds_d"))
    v = verdict(summary)
    cal.to_csv(out / "h31_calibration.tsv", sep="\t", index=False, float_format="%.4f")
    summary.to_csv(out / "h31_summary.tsv", sep="\t", index=False, float_format="%.3f")
    v.to_csv(out / "h31_verdict.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 230, "display.max_colwidth", 110, "display.max_rows", 400):
        print(summary[summary.test.isin(["T1", "T2"])].to_string(index=False)); print(v.to_string(index=False))
    print("wrote", out / "kron_scores.tsv", out / "h31_summary.tsv", out / "h31_verdict.tsv")
    return v


# --- the variant test -------------------------------------------------------------------------------------

def kinds(gb: Path) -> tuple[str, np.ndarray, list]:
    """(the human genome, the kind of each place, the protein genes as `variant_stab.reference` gives them). A place
    in a protein gene is "protein" whatever else overlaps it; then tRNA, rRNA, the control region; else "other"."""
    seq, cds = variant_stab.reference(gb)
    out = np.array(["other"] * len(seq), dtype=object)
    st = stretches(gb, cybrid.GENOMES[cybrid.HUMAN])
    for name, prefix in (("control region", "dloop"), ("rRNA", "rrna"), ("tRNA", "trna"), ("protein", "cds")):
        for k, (_, places, _) in st.items():
            if k[0] == prefix:
                out[places] = name
    return seq, out, cds


def read_variants(clinvar: Path, gnomad: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    """Every single-base variant of H29's two files, whatever it does: harmful (ClinVar, reviewed by the expert panel
    as pathogenic or likely pathogenic) and harmless (gnomAD, passing, homoplasmic in at least AF_HOM_MIN of genomes,
    and pathogenic in no ClinVar record). Columns pos, ref, alt, label, af_hom."""
    counts = {"harmless left out: pathogenic in ClinVar": 0, "harmful without a sequence change given": 0}
    rows, any_pathogenic = {}, set()
    for r in json.loads(clinvar.read_text()):
        if r["type"] != "single nucleotide variant" or not variant_stab.is_pathogenic(r["classification"]):
            continue
        sp = (r["spdi"] or "").split(":")
        if len(sp) != 4 or len(sp[2]) != 1 or len(sp[3]) != 1:
            counts["harmful without a sequence change given"] += r["review"] == variant_stab.PANEL
            continue
        k = (int(sp[1]) + 1, sp[2], sp[3])
        any_pathogenic.add(k)
        if r["review"] == variant_stab.PANEL:
            rows[k] = ("harmful", np.nan)
    g = pd.read_csv(gnomad, sep="\t")
    g = g[(g.filters == "PASS") & (g.ref.str.len() == 1) & (g.alt.str.len() == 1) & (g.AF_hom >= variant_stab.AF_HOM_MIN)]
    for r in g.itertuples():
        k = (int(r.position), r.ref, r.alt)
        if k in any_pathogenic:
            counts["harmless left out: pathogenic in ClinVar"] += 1
        else:
            rows.setdefault(k, ("harmless", float(r.AF_hom)))
    out = pd.DataFrame([{"pos": k[0], "ref": k[1], "alt": k[2], "label": v[0], "af_hom": v[1]} for k, v in rows.items()])
    return out.sort_values(["pos", "alt"]).reset_index(drop=True), counts


def bases_held(data: Path, gb: Path, cls: str = "Mammalia") -> tuple[np.ndarray, list[str]]:
    """([n] how many different bases the build's species of the class hold at each place of the human genome, NaN
    where the build has no alignment of the place; what was not assessed and why). Read from the codon alignments of
    the protein genes and the two rRNA alignments; a gene whose human row is not the record's gene is not assessed."""
    from ..seqio import read_fasta
    seq, _, _ = kinds(gb)
    out, notes = np.full(len(seq), np.nan), ["tRNA genes and the control region: the build has no alignment"]
    sp = pd.read_parquet(data / "species.parquet", columns=["slug", "class"])
    slugs = set(sp.slug[sp["class"] == cls])
    files = {**{("cds", g): f"{g}.codon.fasta.gz" for g in cybrid.ALL_GENES},
             ("rrna", "12S"): "rRNA-12S.nt.fasta.gz", ("rrna", "16S"): "rRNA-16S.nt.fasta.gz"}
    st = stretches(gb, cybrid.GENOMES[cybrid.HUMAN])
    for key, name in files.items():
        aln = read_fasta(data / "alignments" / name)
        human = aln[cybrid.HUMAN].upper()
        gene, places, _ = st[key]
        cols = [c for c, b in enumerate(human) if b != "-"]
        own = "".join(human[c] for c in cols)
        if not gene.startswith(own) and not own.startswith(gene):
            notes.append(f"{key[1]}: the build's human row is not the record's gene")
            continue
        rows = [aln[s].upper() for s in slugs if s in aln]
        if len(rows) < 50:
            notes.append(f"{key[1]}: {len(rows)} species of {cls} in the alignment, fewer than 50")
            continue
        for i, c in enumerate(cols[:len(gene)]):
            n = len({r[c] for r in rows} & set(BASES))
            out[places[i]] = n if np.isnan(out[places[i]]) else min(out[places[i]], n)   # ATP8/ATP6 overlap: the stricter
    return out, notes


def variants(gb: Path, source: Path, data: Path, tables: Path, out: Path, arm: str = "mito") -> pd.DataFrame:
    """kron_variants.tsv (one row per variant: its kind of place, its effect in a protein gene, the score and the
    baseline) and h31_variants_summary.tsv (the rule)."""
    out.mkdir(parents=True, exist_ok=True)
    seq, kind, cds = kinds(gb)
    genome, lp = read_table(tables / f"{cybrid.GENOMES[cybrid.HUMAN]}.{arm}.tsv")
    if genome != seq:
        raise SystemExit("the table is not of the human record's genome")
    v, counts = read_variants(source / "clinvar_mt.json", source / "gnomad.chrM.tsv")
    held, notes = bases_held(data, gb)
    wrong = [f"m.{r.pos}{r.ref}>{r.alt}" for r in v.itertuples() if seq[r.pos - 1] != r.ref]
    v = v[[seq[r.pos - 1] == r.ref for r in v.itertuples()]].copy()
    v["kind"] = [kind[p - 1] for p in v.pos]
    v["effect"] = [";".join(sorted({c[0] for c in variant_stab.consequence(seq, cds, r.pos, r.ref, r.alt)})) for r in v.itertuples()]
    v["kron"] = [drop(lp, r.pos - 1, r.ref, r.alt) for r in v.itertuples()]
    v["n_bases"] = [held[p - 1] for p in v.pos]
    v.to_csv(out / "kron_variants.tsv", sep="\t", index=False, float_format="%.5f")
    rows = [{"what": "variants", "group": "all", "n_harmful": int((v.label == "harmful").sum()), "n_harmless": int((v.label == "harmless").sum()),
             "value": np.nan, "reads": "; ".join([f"{k}: {n}" for k, n in counts.items()] + notes
                                                + [f"reference base not the genome's, left out: {len(wrong)}"])}]
    ok = v[v.kron.notna()].assign(harmful=lambda d: d.label == "harmful")

    def area(d, group):
        h = d.harmful.to_numpy()
        row = {"what": "area, kron", "group": group, "n_harmful": int(h.sum()), "n_harmless": int((~h).sum())}
        if h.sum() < 5 or (~h).sum() < 5:
            rows.append({**row, "value": np.nan, "reads": "not assessed: fewer than 5 of a kind"}); return None
        b = d[d.n_bases.notna()]
        hb = b.harmful.to_numpy()
        rows.append({**row, "value": cybrid.auc(d.kron, d.harmful), "reads": "area under the ROC curve, harmful against harmless"})
        if hb.sum() < 5 or (~hb).sum() < 5:
            rows.append({**row, "what": "difference, kron - baseline", "n_harmful": int(hb.sum()), "n_harmless": int((~hb).sum()),
                         "value": np.nan, "reads": "not assessed: the baseline covers fewer than 5 of a kind"})
            return cybrid.auc(d.kron, d.harmful), None
        a, ab, lo, hi = variant_stab.compare(hb, b.kron.to_numpy(), -b.n_bases.to_numpy())
        rows.append({**row, "what": "area, baseline", "n_harmful": int(hb.sum()), "n_harmless": int((~hb).sum()), "value": ab,
                     "reads": "fewer different bases among the build's mammals is the prediction of harm; where the build has an alignment"})
        rows.append({**row, "what": "difference, kron - baseline", "n_harmful": int(hb.sum()), "n_harmless": int((~hb).sum()), "value": a - ab,
                     "reads": f"on the variants the baseline covers (kron there {a:.3f}); 95% bootstrap interval {lo:.3f} to {hi:.3f}"})
        return cybrid.auc(d.kron, d.harmful), lo

    got = area(ok, "all")
    if got is None:
        verdict_ = "not assessed"
    else:
        verdict_ = "counts" if got[0] >= AUC_MIN and got[1] is not None and got[1] > 0 else "does not count"
    rows.append({"what": "rule", "group": "all", "n_harmful": int(ok.harmful.sum()), "n_harmless": int((~ok.harmful).sum()), "value": np.nan,
                 "reads": f"{verdict_}: an area of at least {AUC_MIN} over all variants and above the baseline with the interval above 0"})
    for k in KINDS:
        area(ok[ok.kind == k], k)
    for e in ("missense", "synonymous"):
        area(ok[ok.effect == e], f"protein, {e}")
    res = pd.DataFrame(rows)
    res.to_csv(out / "h31_variants_summary.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 230, "display.max_colwidth", 120, "display.max_rows", 200):
        print(res.to_string(index=False))
    print("wrote", out / "kron_variants.tsv", out / "h31_variants_summary.tsv")
    return res
