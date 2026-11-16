"""R2 (H21-H23) and L3 (H24): one pair per origin of a trait, each sister against one outgroup.

Per pair (`clades.origin_pairs`): dN and dS of the trait sister and of the other sister to the outgroup, over the 12
heavy-strand genes concatenated (Nei-Gojobori with the Jukes-Cantor correction, in the species' mitochondrial code, as
H7 counts them). Both sisters have had the same time since they split from the outgroup's line, so a difference of
log10 dS between them is also a difference of log10 dS per million years, and needs no dating.

    d_omega = log10(dN/dS of the trait sister) - log10(dN/dS of the other)
    d_dS    = log10 dS(trait, outgroup) - log10 dS(other, outgroup);  d_dN likewise

A pair is unusable when any of the four distances is missing or 0, and saturated when either dS to the outgroup is
above 1; saturated pairs are left out of the primary test and reported with it as a sensitivity.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from scipy import stats

from ..genes import MT_HEAVY_STRAND
from ..longevity.rates import SATURATED_DS, distance
from ..seqio import read_fasta
from . import clades

MIN_PAIRS = 8
BOOT = 10_000


class Distances:
    """dS and dN between two species over the 12 heavy-strand genes concatenated, computed once per pair."""
    def __init__(self, data: Path, species: pd.DataFrame, genes=MT_HEAVY_STRAND):
        self.genes = list(genes)
        self.aln = {g: read_fasta(data / "alignments" / f"{g}.codon.fasta.gz") for g in self.genes}
        self.aligned = set.intersection(*[set(a) for a in self.aln.values()])
        self.table = species.set_index("slug").transl_table
        self.cache: dict[tuple[str, str], dict] = {}

    def code(self, a: str, b: str) -> int:
        tables = {int(self.table[s]) for s in (a, b) if not pd.isna(self.table.get(s))}
        t = tables.pop() if len(tables) == 1 else 2
        return t if t in CodonTable.unambiguous_dna_by_id else 2

    def __call__(self, a: str, b: str) -> dict:
        key = (a, b) if a <= b else (b, a)
        if key not in self.cache:
            sa, sb = ("".join(self.aln[g][s] for g in self.genes) for s in key)
            self.cache[key] = distance(sa, sb, self.code(*key))
        return self.cache[key]


def pair_rates(rows: list[dict], dist: Distances) -> pd.DataFrame:
    out = []
    for r in rows:
        if not r["taken"]:
            out.append(r); continue
        t, o, g = r["sp_trait"], r["sp_other"], r["outgroup"]
        dt, do = dist(t, g), dist(o, g)
        vals = {"dS_trait": dt["dS"], "dN_trait": dt["dN"], "dS_other": do["dS"], "dN_other": do["dN"]}
        ok = all(np.isfinite(v) and v > 0 for v in vals.values())
        row = {**r, **vals, "usable": ok, "saturated": bool(ok and max(vals["dS_trait"], vals["dS_other"]) > SATURATED_DS)}
        if ok:
            row["d_dS"] = np.log10(vals["dS_trait"]) - np.log10(vals["dS_other"])
            row["d_dN"] = np.log10(vals["dN_trait"]) - np.log10(vals["dN_other"])
            row["d_omega"] = row["d_dN"] - row["d_dS"]
        out.append(row)
    return pd.DataFrame(out)


def median_ci(x: np.ndarray, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    b = np.median(rng.choice(x, size=(BOOT, len(x)), replace=True), axis=1)
    return float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))


def test(p: pd.DataFrame, col: str, with_saturated: bool = False) -> dict:
    """Wilcoxon signed-rank (two-sided) of `col` over the usable pairs, with the sign test, the median and its
    bootstrap 95% interval."""
    d = p[p.get("usable", pd.Series(False, index=p.index)).fillna(False).astype(bool)]
    if not with_saturated:
        d = d[~d.saturated.astype(bool)]
    x = d[col].to_numpy(dtype=float)
    base = {"measure": col, "with_saturated": with_saturated, "n_pairs": len(x)}
    if len(x) < MIN_PAIRS:
        return {**base, "status": f"not assessed: {len(x)} pairs, fewer than {MIN_PAIRS}"}
    w = stats.wilcoxon(x)
    pos = int((x > 0).sum()); nz = int((x != 0).sum())
    lo, hi = median_ci(x)
    return {**base, "status": "assessed", "median": float(np.median(x)), "ci_low": lo, "ci_high": hi,
            "p": float(w.pvalue), "positive": pos, "sign_p": float(stats.binomtest(pos, nz).pvalue) if nz else np.nan}


def size_adjusted(p: pd.DataFrame, size: pd.Series, col: str = "d_omega") -> dict:
    """The intercept of `col` regressed on the pair's difference in size (trait sister minus other), OLS."""
    d = p[p.get("usable", False).fillna(False).astype(bool) & ~p.saturated.fillna(False).astype(bool)].copy()
    d["d_size"] = d.sp_trait.map(size) - d.sp_other.map(size)
    d = d.dropna(subset=["d_size", col])
    base = {"measure": f"{col} given size", "n_pairs": len(d)}
    if len(d) < MIN_PAIRS:
        return {**base, "status": f"not assessed: {len(d)} pairs with both sizes"}
    X = np.column_stack([np.ones(len(d)), d.d_size.to_numpy(float)])
    y = d[col].to_numpy(float)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b; df = len(y) - 2
    se = np.sqrt(np.diag(np.linalg.inv(X.T @ X)) * (r @ r) / df)
    t = b[0] / se[0]
    return {**base, "status": "assessed", "intercept": b[0], "se": se[0], "p": 2 * stats.t.sf(abs(t), df),
            "size_slope": b[1]}


def run(tree_file: Path, traits: pd.DataFrame, trait: str, usable: set[str], dist: Distances,
        scope: pd.Series | None = None) -> pd.DataFrame:
    """Pairs of one trait on the dated tree pruned to the species with a known state that are `usable`."""
    t = traits[traits.trait.eq(trait) & traits.state.notna()]
    if scope is not None:
        t = t[t.slug.isin(set(scope))]
    state = {s: int(v) for s, v in zip(t.slug, t.state) if s in usable}
    root = clades.prune(clades.read(tree_file), set(state))
    if root is None:
        return pd.DataFrame()
    state = {s: state[s] for s in clades.tips(root)}
    return pair_rates(clades.origin_pairs(root, state), dist).assign(trait=trait)
