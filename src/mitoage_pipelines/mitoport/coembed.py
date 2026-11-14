"""Do a nuclear protein's ESM-2 embedding distances between species follow those of the mitochondrial
subunits of a complex, beyond what the dated species tree gives?

For a nuclear protein g and a mitochondrial set M (mean-pooled ESM-2 650M embeddings, cosine distance
between species; for M the distances of its proteins are averaged after scaling each to mean 1): partial
correlation over species pairs of d_g and d_M given the patristic distance on the dated tree and its
logarithm (ranks, then residuals). Species pairs are not independent, so no p-value is given per protein:
groups of proteins are compared with the cytosolic controls (Mann-Whitney over proteins).
"""
from pathlib import Path

import numpy as np
import polars as pl
from Bio import Phylo
from scipy import stats

from ..genes import GROUPS

MT = {"IV": ["COX1", "COX2", "COX3"], "III": ["CYTB"], "I": "ND1 ND2 ND3 ND4 ND4L ND5 ND6".split(), "V": ["ATP6", "ATP8"]}
MIN_SPECIES = 60


def load(folder: Path, prefix: str) -> dict[str, dict[str, np.ndarray]]:
    """{gene: {species slug: embedding}} from <prefix>.index.parquet (ids GENE|slug) and <prefix>.mean.npy."""
    idx = pl.read_parquet(folder / f"{prefix}.index.parquet")
    X = np.load(folder / f"{prefix}.mean.npy").astype(np.float64)
    out: dict[str, dict[str, np.ndarray]] = {}
    for i, pid in enumerate(idx["id"]):
        gene, slug = pid.split("|")
        out.setdefault(gene, {})[slug] = X[i]
    return out


def cosine_distance(A: np.ndarray) -> np.ndarray:
    A = A / np.linalg.norm(A, axis=1, keepdims=True)
    return 1 - A @ A.T


def residual(y: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    Z = np.column_stack([np.ones(len(y))] + covariates)
    b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return y - Z @ b


def patristic(tree, tips: dict, depth: dict, species: list[str]) -> np.ndarray:
    n = len(species); D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a = tree.common_ancestor(tips[species[i]], tips[species[j]])
            D[i, j] = depth[tips[species[i]]] + depth[tips[species[j]]] - 2 * depth[a]
    return D


def main(embeddings: Path, tree_file: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    nuc, mt = load(embeddings, "nuclear-esm2_650m"), load(embeddings, "mt13-esm2_650m")
    groups = dict(GROUPS)
    known = {g for v in groups.values() for g in v}
    groups["cytosolic control"] = sorted(g for g in nuc if g not in known)
    tree = Phylo.read(tree_file, "newick")
    tips = {t.name: t for t in tree.get_terminals()}; depth = tree.depths()
    cache: dict[tuple, np.ndarray] = {}
    rows = []
    for grp, genes in groups.items():
        for g in genes:
            if g not in nuc: continue
            sp = sorted(s for s in nuc[g] if s in tips and all(s in mt[m] for ms in MT.values() for m in ms))
            if len(sp) < MIN_SPECIES: continue
            iu = np.triu_indices(len(sp), 1)
            key = tuple(sp)
            if key not in cache: cache[key] = patristic(tree, tips, depth, sp)
            T = cache[key][iu]
            rg = residual(stats.rankdata(cosine_distance(np.array([nuc[g][s] for s in sp]))[iu]), [T, np.log1p(T)])
            r = {"group": grp, "gene": g, "species": len(sp)}
            for cx, ms in MT.items():
                ds = [cosine_distance(np.array([mt[m][s] for s in sp]))[iu] for m in ms]
                dm = np.mean([d / d.mean() for d in ds], axis=0)
                r[cx] = float(np.corrcoef(rg, residual(stats.rankdata(dm), [T, np.log1p(T)]))[0, 1])
            rows.append(r)
    res = pl.DataFrame(rows)
    res.write_csv(out / "coembed.tsv", separator="\t", float_precision=3)
    ctrl = res.filter(pl.col("group") == "cytosolic control")
    print("median partial correlation with each mitochondrial complex; Mann-Whitney p against the cytosolic controls")
    for grp in groups:
        d = res.filter(pl.col("group") == grp)
        if d.height == 0: continue
        line = f"{grp:28s} n={d.height:2d}"
        for cx in MT:
            p = stats.mannwhitneyu(d[cx], ctrl[cx])[1] if grp != "cytosolic control" and d.height > 1 else float("nan")
            line += f"  {cx} {d[cx].median():+.3f} (p {p:.2g})"
        print(line)
