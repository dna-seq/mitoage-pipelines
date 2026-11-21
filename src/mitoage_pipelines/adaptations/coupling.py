"""C, the coevolution side of H21's design (H21, H23; H23's positive control P).

Species: the class's species that are tips of the dated tree and have the 13 mitochondrial proteins and every protein
of H19's vertebrate alignments of the nuclear subunits of complexes I, III, IV and V and of the 30 cytosolic controls.
If fewer than 5 trait-bearing species remain, the proteins kept are those present in at least 90% of the class's
species with the 13, and the species those that have all of them; `prepare` says which rule it used. When no species
is left (the nuclear alignments cover fewer than 90% of the class), `run` answers "not assessed" and no tree is fitted.

`prepare` writes, per class, one alignment per protein restricted to those species (columns left all-gap dropped),
the dated topology pruned to them, genes.txt and species.txt; `mitoport branch-lengths` (IQ-TREE, H5's models) then
gives every protein's branch lengths on that topology.

Per branch b (`erc.branches`: keyed by the smaller side of its bipartition): each protein's length scaled by its mean
over the branches; the nuclear background, the mean over the cytosolic controls of their scaled lengths (a control's own
left out of its background); a nuclear protein's relative rate r = scaled / background; the mitochondrial rate m = the
mean scaled length of the 13 / background. A length of 0 has no logarithm: that branch is left out of that protein's
medians and counted.
- C1: per protein, s = median over foreground of log r - median over the other branches; the subunits of the four
  complexes against the cytosolic controls, Mann-Whitney, two-sided.
- C2: median over foreground of log m - median over the other branches, against matched control clade sets
  (`clades.control_draws`), two-sided empirical p.
- P (H23): C1 with the bats' stem branch alone as foreground, in mammals.
Foreground: the stem and every inner branch of each trait clade, as clades of the topology pruned to the species set.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ..genes import COMPLEX_I_NUCLEAR, COMPLEX_III_NUCLEAR, COMPLEX_IV_NUCLEAR, COMPLEX_V_NUCLEAR, CYTOSOLIC_CONTROLS, MT_GENES
from ..mitoport.erc import branches
from ..seqio import read_fasta
from . import clades

SUBUNITS = COMPLEX_I_NUCLEAR + COMPLEX_III_NUCLEAR + COMPLEX_IV_NUCLEAR + COMPLEX_V_NUCLEAR
MIN_TRAIT_SPECIES = 5
MIN_TRAIT_CLADES = 3
RELAXED_SHARE = 0.9
N_NULL = 1000


def mt_name(g: str) -> str:
    return f"MT-{g}"


def prepare(data: Path, nucaln: Path, tree_file: Path, species: pd.DataFrame, cls: str, trait_species: set[str],
            out: Path) -> dict:
    tree = clades.read(tree_file)
    tipset = set(clades.tips(tree))
    mt = {g: read_fasta(data / "alignments" / f"{g}.aa.fasta.gz") for g in MT_GENES}
    nuc = {g: read_fasta(nucaln / f"{g}.afa") for g in SUBUNITS + CYTOSOLIC_CONTROLS if (nucaln / f"{g}.afa").is_file()}
    missing = [g for g in SUBUNITS + CYTOSOLIC_CONTROLS if g not in nuc]
    base = set(species[species["class"].eq(cls) & species.record_status.ne("retired")].slug) & tipset
    base &= set.intersection(*[{s for s, v in a.items() if v.strip("-")} for a in mt.values()])
    has = {g: {s for s, v in a.items() if v.strip("-")} & base for g, a in nuc.items()}
    genes = sorted(nuc)
    sp = base.intersection(*has.values()) if has else set()
    rule = "every protein"
    if len(sp & trait_species) < MIN_TRAIT_SPECIES:
        genes = sorted(g for g in nuc if len(has[g]) >= RELAXED_SHARE * len(base))
        sp = base.intersection(*[has[g] for g in genes]) if genes else set()
        rule = f"proteins in at least {RELAXED_SHARE:.0%} of the class's species"
    out.mkdir(parents=True, exist_ok=True)
    order = sorted(sp)
    allg = [mt_name(g) for g in MT_GENES] + genes
    for g in allg:
        a = mt[g[3:]] if g.startswith("MT-") else nuc[g]
        rows = [a[s].upper() for s in order]
        keep = [i for i in range(len(rows[0])) if any(r[i] != "-" for r in rows)] if rows else []
        (out / f"{g}.faa").write_text("".join(f">{s}\n{''.join(r[i] for i in keep)}\n" for s, r in zip(order, rows)))
    root = clades.prune(tree, sp)
    (out / "topology.nwk").write_text(clades.newick(root, raw_names=True) if root else "")
    (out / "genes.txt").write_text("\n".join(allg) + "\n")
    (out / "species.txt").write_text("\n".join(order) + "\n")
    info = {"class": cls, "rule": rule, "class_species": len(base),
            "with_a_nuclear_protein": len(set().union(*has.values())) if has else 0,
            "species": len(order), "trait_species": len(sp & trait_species),
            "nuclear_proteins": len(genes), "subunits": sum(g in SUBUNITS for g in genes),
            "controls": sum(g in CYTOSOLIC_CONTROLS for g in genes), "no_alignment": " ".join(missing)}
    pd.DataFrame([info]).to_csv(out / "prepare.tsv", sep="\t", index=False)
    return info


def relative_rates(bl: dict[str, dict[frozenset, float]], genes: list[str]) -> tuple[pd.DataFrame, pd.Series, dict]:
    """-> (log r per branch x nuclear protein, log m per branch, branches left out per protein for a length of 0)."""
    keys = sorted(set.intersection(*[set(bl[g]) for g in genes]), key=lambda k: (len(k), sorted(k)))
    L = pd.DataFrame({g: [bl[g][k] for k in keys] for g in genes}, index=pd.Index(range(len(keys))))
    scaled = L / L.mean()
    ctl = [g for g in genes if g in CYTOSOLIC_CONTROLS]
    mts = [g for g in genes if g.startswith("MT-")]
    nucl = [g for g in genes if not g.startswith("MT-")]
    tot = scaled[ctl].sum(axis=1)
    logr, zeros = {}, {}
    for g in nucl:
        bg = (tot - scaled[g]) / (len(ctl) - 1) if g in ctl else tot / len(ctl)
        r = scaled[g] / bg
        zeros[g] = int((r <= 0).sum())
        logr[g] = np.log(r.where(r > 0))
    m = scaled[mts].mean(axis=1) / (tot / len(ctl))
    out = pd.DataFrame(logr); out.attrs["keys"] = keys
    return out, np.log(m.where(m > 0)), zeros


def foreground(keys: list[frozenset], full: clades.Node, cl: list[clades.Node], known: set[str],
               stems_only: bool = False) -> np.ndarray:
    """Branches of the species-set topology `full` whose known tips (in `known`) are not empty and all lie in one
    clade of `cl` (a tip of unknown state is ignored, so it neither joins nor breaks a clade); with `stems_only`, the
    topmost such branch of each clade."""
    alltips = frozenset(clades.tips(full))
    norm = lambda z: min(z, alltips - z, key=lambda x: (len(x), sorted(x)))
    sets = [frozenset(clades.tips(c)) for c in cl]
    inside = {}
    for n in clades.walk(full):
        k = frozenset(clades.tips(n)) & known
        inside[id(n)] = bool(k) and any(k <= cs for cs in sets)
    want = set()
    for n in clades.walk(full):
        if n.parent is None or not inside[id(n)]:
            continue
        if stems_only and inside[id(n.parent)]:
            continue
        want.add(norm(frozenset(clades.tips(n))))
    return np.array([k in want for k in keys])


def c1(logr: pd.DataFrame, fg: np.ndarray, readout: str) -> dict:
    s = {g: logr[g][fg].median() - logr[g][~fg].median() for g in logr}
    sub = pd.Series({g: v for g, v in s.items() if g in SUBUNITS}).dropna()
    ctl = pd.Series({g: v for g, v in s.items() if g in CYTOSOLIC_CONTROLS}).dropna()
    base = {"readout": readout, "foreground_branches": int(fg.sum()), "n_subunits": len(sub), "n_controls": len(ctl)}
    if len(sub) < 3 or len(ctl) < 3 or fg.sum() == 0:
        return {**base, "status": "not assessed: too few proteins or no foreground branch"}
    u = stats.mannwhitneyu(sub, ctl, alternative="two-sided")
    return {**base, "status": "assessed", "median_subunits": sub.median(), "median_controls": ctl.median(),
            "estimate": sub.median() - ctl.median(), "p": float(u.pvalue)}


def c2(logm: pd.Series, fg: np.ndarray, null_fgs: list[np.ndarray]) -> dict:
    logm = logm.where(np.isfinite(logm))
    stat = lambda f: logm[f].median() - logm[~f].median()
    obs = stat(fg)
    null = np.array([stat(f) for f in null_fgs if f.any()])
    null = null[np.isfinite(null)]
    base = {"readout": "C2", "foreground_branches": int(fg.sum()), "null_draws": len(null),
            "invalid_null_draws": len(null_fgs) - len(null)}
    if not np.isfinite(obs):
        return {**base, "status": "not assessed: no finite foreground/background contrast"}
    if len(null) < 100 or not fg.any():
        return {**base, "status": "not assessed: too few control draws or no foreground branch"}
    mid = np.median(null)
    p = (1 + (np.abs(null - mid) >= abs(obs - mid)).sum()) / (1 + len(null))
    return {**base, "status": "assessed", "estimate": obs, "null_median": mid, "p": float(p)}


def run(work: Path, bl_dir: Path, traits: pd.DataFrame, trait: str, positive_control: bool = False,
        n_null: int = N_NULL) -> list[dict]:
    genes = (work / "genes.txt").read_text().split()
    sp = set((work / "species.txt").read_text().split())
    if not sp:      # `prepare` kept no species: said as such, never run on empty alignments
        return [{"readout": r, "trait": trait, "status": "not assessed: no species has every protein the rule keeps"}
                for r in (("P",) if positive_control else ("C1", "C2"))]
    bl = {}
    for g in genes:
        f = bl_dir / f"{g}.treefile"
        if not f.is_file():
            return [{"readout": r, "trait": trait, "status": f"not assessed: no branch lengths for {g}"} for r in ("C1", "C2")]
        bl[g] = branches(f)
    logr, logm, zeros = relative_rates(bl, genes)
    keys = logr.attrs["keys"]
    t = traits[traits.trait.eq(trait) & traits.state.notna() & traits.slug.isin(sp)]
    state = {s: int(v) for s, v in zip(t.slug, t.state)}
    root = clades.prune(clades.read(work / "topology.nwk"), set(state))
    state = {s: state[s] for s in clades.tips(root)}
    tc = clades.trait_clades(root, state)
    out = []
    n_tr = sum(v == 1 for v in state.values())
    if positive_control:
        full = clades.prune(clades.read(work / "topology.nwk"), sp)
        r = c1(logr, foreground(keys, full, tc, set(state), stems_only=True), "P")
        return [{**r, "trait": trait, "trait_clades": len(tc), "zero_length_branches": sum(zeros.values())}]
    if len(tc) < MIN_TRAIT_CLADES or n_tr < MIN_TRAIT_SPECIES:
        return [{"readout": r, "trait": trait, "trait_clades": len(tc), "trait_species": n_tr,
                 "status": f"not assessed: {len(tc)} trait clades and {n_tr} trait species in the species set"}
                for r in ("C1", "C2")]
    full = clades.prune(clades.read(work / "topology.nwk"), sp)
    fg = foreground(keys, full, tc, set(state))
    draws, short = clades.control_draws(root, state, tc, n_null)
    null = [foreground(keys, full, d, set(state)) for d in draws]
    for r in (c1(logr, fg, "C1"), c2(logm, fg, null)):
        out.append({**r, "trait": trait, "trait_clades": len(tc), "trait_species": n_tr,
                    "zero_length_branches": sum(zeros.values()), "draws_short": int(sum(s > 0 for s in short))})
    return out
