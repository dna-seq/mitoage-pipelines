"""R1 (RELAX) and R4 (aBSREL) of H21's design, for H21-H23 and H25's B4.

Inputs, per trait and draw (draw 0 = the trait clades; draws 1..20 = matched control clades, `clades.control_draws`):
- tips: every species of the draw's clades, and for each clade the three state-0 species nearest its root, outside
  every clade of the draw;
- the dated tree pruned to them; for RELAX every branch inside a clade and its stem is labelled {Test}, every other
  {Reference}; for aBSREL only the stems are {Test} (tested), the rest unlabelled;
- the 12 heavy-strand genes concatenated as codons: per gene, codon columns filled in at least half of the tips are
  kept, and a stop, partial or ambiguous codon is written as a gap. Every tip must use the vertebrate mitochondrial
  code (table 2), which HyPhy is told; a tip with another code is left out and counted.

HyPhy 2.5.104: `relax --models Minimal` (the test alone) and `absrel`. A JSON without the fields read here is refused,
never read as a null result.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Data import CodonTable

from .. import external
from ..genes import MT_HEAVY_STRAND
from ..seqio import read_fasta
from . import clades

N_REFERENCE = 3
N_DRAWS = 20
CODE = "Vertebrate-mtDNA"


class Refused(RuntimeError):
    pass


def draw_tips(root: clades.Node, state: dict[str, int], cl: list[clades.Node]) -> set[str]:
    inside = set().union(*[set(clades.tips(c)) for c in cl]) if cl else set()
    keep = set(inside)
    for c in cl:
        n, got = c, []
        while n.parent is not None and len(got) < N_REFERENCE:
            a = n.parent
            cands = sorted((clades.depth_to(t, a), t.name) for t in clades.walk(a)
                           if t.is_tip and t.name not in inside and state[t.name] == 0 and t.name not in got)
            got += [s for _, s in cands[:N_REFERENCE - len(got)]]
            n = a
        keep |= set(got)
    return keep


def labelled(root: clades.Node, cl: list[clades.Node], keep: set[str], stems_only: bool) -> tuple[clades.Node, str, int]:
    """The tree pruned to `keep`, as newick with the {Test} labels; -> (tree, newick, number of Test branches)."""
    sets = [frozenset(clades.tips(c)) for c in cl]
    t = clades.prune(root, keep)
    test = set()
    for n in clades.walk(t):
        s = frozenset(clades.tips(n))
        for cs in sets:
            if s == cs or (not stems_only and s <= cs):
                test.add(id(n))
    if stems_only:
        nwk = clades.newick(t, None)
        def rec(n):
            lab = "{Test}" if id(n) in test else ""
            if n.is_tip:
                return f"{n.name.replace('-', '_')}{lab}:{n.length:.6g}"
            return "(" + ",".join(rec(c) for c in n.children) + f"){lab}:{n.length:.6g}"
        nwk = rec(t) + ";\n"
    else:
        nwk = clades.newick(t, test)
    return t, nwk, len(test)


def write_codons(data: Path, tips: list[str], out: Path) -> dict:
    stops = set(CodonTable.unambiguous_dna_by_id[2].stop_codons)
    seqs = {s: [] for s in tips}; kept = {}; masked = 0
    for g in MT_HEAVY_STRAND:
        aln = read_fasta(data / "alignments" / f"{g}.codon.fasta.gz")
        rows = {s: aln[s].upper() for s in tips}
        L = len(next(iter(rows.values())))
        if L % 3:
            raise Refused(f"{g}: codon alignment length {L} is not a multiple of 3")
        keep = [i for i in range(0, L, 3) if sum(r[i:i + 3] != "---" for r in rows.values()) >= 0.5 * len(rows)]
        kept[g] = len(keep)
        for s, r in rows.items():
            for i in keep:
                c = r[i:i + 3]
                if c in stops or ("-" in c and c != "---") or any(x not in "ACGT-" for x in c):
                    masked += c != "---"; c = "---"
                seqs[s].append(c)
    out.write_text("".join(f">{s.replace('-', '_')}\n{''.join(v)}\n" for s, v in seqs.items()))
    return {"codons": sum(kept.values()), "masked": masked, **{f"codons_{g}": n for g, n in kept.items()}}


def prepare(data: Path, tree_file: Path, traits: pd.DataFrame, trait: str, usable: set[str], codes: pd.Series,
            out: Path, scope: set[str] | None = None, n_draws: int = N_DRAWS) -> pd.DataFrame:
    """Writes <out>/<draw>/{tree.nwk, absrel.nwk, concat.fna, tips.txt} for draw 0 (the trait) and the controls."""
    t = traits[traits.trait.eq(trait) & traits.state.notna()]
    if scope is not None:
        t = t[t.slug.isin(scope)]
    wrong = sorted(s for s in t.slug if s in usable and codes.get(s) != 2)
    state = {s: int(v) for s, v in zip(t.slug, t.state) if s in usable and codes.get(s) == 2}
    root = clades.prune(clades.read(tree_file), set(state))
    state = {s: state[s] for s in clades.tips(root)}
    tc = clades.trait_clades(root, state)
    draws, short = clades.control_draws(root, state, tc, n_draws)
    rows = []
    for k, cl in enumerate([tc] + draws):
        d = out / f"draw{k:02d}"; d.mkdir(parents=True, exist_ok=True)
        keep = draw_tips(root, state, cl)
        _, nwk, n_test = labelled(root, cl, keep, stems_only=False)
        _, anwk, n_stem = labelled(root, cl, keep, stems_only=True)
        (d / "tree.nwk").write_text(nwk); (d / "absrel.nwk").write_text(anwk)
        (d / "tips.txt").write_text("\n".join(sorted(keep)) + "\n")
        info = write_codons(data, sorted(keep), d / "concat.fna")
        rows.append({"trait": trait, "draw": k, "clades": len(cl), "clades_left_out": 0 if k == 0 else short[k - 1],
                     "tips": len(keep), "test_branches": n_test, "stems": n_stem, "other_code_left_out": len(wrong), **info})
    res = pd.DataFrame(rows)
    res.to_csv(out / "draws.tsv", sep="\t", index=False)
    return res


def run_relax(d: Path, cpus: int = 1) -> Path:
    out = d / "relax.json"
    external.run([external.find("hyphy"), "relax", f"CPU={cpus}", "--alignment", d / "concat.fna", "--tree",
                  d / "tree.nwk", "--code", CODE, "--test", "Test", "--reference", "Reference", "--models", "Minimal",
                  "--output", out])
    return out


def run_absrel(d: Path, cpus: int = 1) -> Path:
    out = d / "absrel.json"
    external.run([external.find("hyphy"), "absrel", f"CPU={cpus}", "--alignment", d / "concat.fna", "--tree",
                  d / "absrel.nwk", "--code", CODE, "--branches", "Test", "--output", out])
    return out


def read_relax(path: Path) -> dict:
    r = json.loads(Path(path).read_text())
    try:
        t = r["test results"]
        return {"LRT": float(t["LRT"]), "p": float(t["p-value"]), "K": float(t["relaxation or intensification parameter"])}
    except (KeyError, TypeError, ValueError) as e:
        raise Refused(f"{path}: not a RELAX result with LRT, p-value and K ({e!r})") from e


def read_absrel(path: Path) -> pd.DataFrame:
    """One row per tested branch: its corrected and uncorrected p (aBSREL's own Holm-Bonferroni over tested branches)."""
    r = json.loads(Path(path).read_text())
    try:
        tested = r["tested"]["0"]
        attrs = r["branch attributes"]["0"]
        rows = [{"branch": b, "p_corrected": float(attrs[b]["Corrected P-value"]),
                 "p_uncorrected": float(attrs[b]["Uncorrected P-value"])}
                for b, v in tested.items() if v == "test"]
    except (KeyError, TypeError, ValueError) as e:
        raise Refused(f"{path}: not an aBSREL result with tested branches and their p-values ({e!r})") from e
    if not rows:
        raise Refused(f"{path}: aBSREL tested no branch")
    return pd.DataFrame(rows)


def summarise(out: Path) -> pd.DataFrame:
    """relax.tsv: K, LRT, p of every draw found; the trait's LRT against the controls'."""
    rows = []
    for d in sorted(out.glob("draw*")):
        f = d / "relax.json"
        rows.append({"draw": int(d.name[4:]), **(read_relax(f) if f.is_file() else {"status": "missing"})})
    res = pd.DataFrame(rows)
    res.to_csv(out / "relax.tsv", sep="\t", index=False)
    return res


def r1(res: pd.DataFrame, n_draws: int = N_DRAWS) -> dict:
    """R1's numbers: the trait's K, LRT and p, and how many of the control draws' LRTs it exceeds."""
    if res.empty or "LRT" not in res or res.loc[res.draw.eq(0), "LRT"].isna().all():
        return {"status": "not assessed: the trait's RELAX result is missing"}
    tr = res[res.draw.eq(0)].iloc[0]
    ctl = res[res.draw.gt(0)].LRT.dropna()
    out = {"K": tr.K, "LRT": tr.LRT, "p": tr.p, "controls_run": int(len(ctl)),
           "controls_below": int((ctl < tr.LRT).sum())}
    if len(ctl) < n_draws:
        return {**out, "status": f"not assessed: {len(ctl)} of {n_draws} control draws have a result"}
    return {**out, "status": "assessed"}
