"""The column test with a tree: phylogenetic regression of log10 lifespan on residue state, per column.

Covariance: shared root-to-ancestor path on the tree (Brownian), scaled by Pagel's lambda, which is
estimated once on the model without the site (lifespan ~ mass + fourfold composition). The tree is the
build's dated species tree; a tree built from the genes (with its outgroup) is accepted only to show how
much the result depends on the tree.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import linalg, optimize, stats
from statsmodels.stats.multitest import multipletests

from ..genes import MT_GENES
from .alignments import protein_alignment


def shared_path_covariance(clade, tips: dict[str, int]) -> np.ndarray:
    """C[i, j] = length of the path from `clade` to the last common ancestor of tips i and j
    (on the diagonal: the tip's distance from `clade`). `tips` maps each tip name under `clade` to its row."""
    C = np.zeros((len(tips), len(tips)))
    def walk(cl, depth):
        for ch in cl.clades:
            d = depth + (ch.branch_length or 0.0)
            ids = [tips[t.name] for t in ch.get_terminals()]
            C[np.ix_(ids, ids)] = d          # top-down: a deeper clade overwrites its own block
            if not ch.is_terminal(): walk(ch, d)
    walk(clade, 0.0)
    return C


def V(lam, C):
    v = lam * C; np.fill_diagonal(v, np.diag(C)); return v


def nll(lam, X, y, C):
    L = linalg.cholesky(V(lam, C), lower=True)
    Xw, yw = linalg.solve_triangular(L, X, lower=True), linalg.solve_triangular(L, y, lower=True)
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw @ b; n = len(y)
    return 0.5 * n * np.log(r @ r / n) + np.log(np.diag(L)).sum()


def gls(X, y, Cs, lam):
    L = linalg.cholesky(V(lam, Cs), lower=True)
    Xw, yw = linalg.solve_triangular(L, X, lower=True), linalg.solve_triangular(L, y, lower=True)
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw @ b; df = len(y) - X.shape[1]
    cov = np.linalg.pinv(Xw.T @ Xw) * (r @ r / df)
    return b, np.sqrt(np.diag(cov)), df


def main(data: Path, composition: Path, fourfold: Path, tree_file: Path, outgroup_file: Path | None, out: Path,
         cls: str = "Mammalia", alt_dir: Path | None = None, alt_tool: str | None = None,
         quality_file: Path | None = None) -> None:
    """`out` is the table to write. `alt_dir`/`alt_tool`: read <alt_dir>/<GENE>.<alt_tool>.afa instead of the build's alignments.
    `outgroup_file` None: `tree_file` is the dated species tree, rooted, and the class is taken as the tips under
    the common ancestor of its species. `quality_file`: keep only the species listed in it."""
    out.parent.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(composition, sep="\t").merge(pd.read_csv(fourfold, sep="\t"), on="slug")
    sp = sp[sp["class"] == cls].dropna(subset=["lm"]).reset_index(drop=True)
    if quality_file is not None:
        okq = set(quality_file.read_text().split()); n0 = len(sp)
        sp = sp[sp.slug.isin(okq)].reset_index(drop=True)
        print(f"lifespan quality high or acceptable: {len(sp)} of {n0}")
    tree = Phylo.read(tree_file, "newick")
    if outgroup_file is None:
        present = {t.name for t in tree.get_terminals()}
        sp = sp[sp.slug.isin(present)].reset_index(drop=True)
        wanted = set(sp.slug)
        anc = tree.common_ancestor([t for t in tree.get_terminals() if t.name in wanted])
        extra = sum(t.name not in wanted for t in anc.get_terminals())
        print(f"dated tree: {len(sp)} species of {cls} in it; {extra} other tips under their common ancestor (ignored)")
    else:
        outg = [l.strip() for l in outgroup_file.read_text().split()]
        tree.root_with_outgroup(*[t for t in tree.get_terminals() if t.name in outg][:1])
        anc = tree.common_ancestor([t for t in tree.get_terminals() if t.name not in outg])
        mono = all(t.name not in outg for t in anc.get_terminals())
        print(f"{cls} monophyletic against the outgroup:", mono)
        if not mono: raise SystemExit(f"outgroup falls inside {cls}: the rooting is wrong, stopping")
    tips = {t.name: i for i, t in enumerate(anc.get_terminals())}
    sp = sp[sp.slug.isin(tips)].reset_index(drop=True)
    idx = np.array([tips[s] for s in sp.slug])
    C = shared_path_covariance(anc, tips)
    C = C[np.ix_(idx, idx)]
    print(f"species: {len(sp)}; root-to-tip height median {np.median(np.diag(C)):.3f}, min {np.diag(C).min():.3f}, max {np.diag(C).max():.3f}")
    Y = sp.ll.to_numpy(); COV = np.column_stack([np.ones(len(sp)), sp[["lm", "A4", "C4", "T4"]].to_numpy()])

    opt = optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(COV, Y, C), method="bounded")
    lam = float(opt.x); print(f"Pagel's lambda of lifespan given mass and composition: {lam:.3f}")
    b, se, df = gls(COV, Y, C, lam)
    print("null model: mass slope %.3f (t %.1f)" % (b[1], b[1] / se[1]))

    fam = sp.family.to_numpy(); rows = []
    for g in MT_GENES:
        aln = protein_alignment(data, g, alt_dir, alt_tool)
        M = np.array([list(aln[s]) for s in sp.slug]); hrow = aln["homo-sapiens"]; hpos = np.cumsum([c != "-" for c in hrow])
        for col in range(M.shape[1]):
            c = M[:, col]; ok = c != "-"
            if ok.mean() < 0.8: continue
            Cs = C[np.ix_(ok, ok)] if not ok.all() else C
            for aa in np.unique(c[ok]):
                x = (c[ok] == aa).astype(float)
                if x.sum() < 15 or (1 - x).sum() < 15: continue
                d = pd.DataFrame({"f": fam[ok], "x": x}).groupby("f").x.agg(["sum", "size"])
                poly = int(((d["sum"] > 0) & (d["sum"] < d["size"])).sum())
                if poly < 5: continue
                b, se, df = gls(np.column_stack([COV[ok], x]), Y[ok], Cs, lam)
                t = b[-1] / se[-1]
                rows.append(dict(gene=g, column=col + 1, human_pos=int(hpos[col]) if hrow[col] != "-" else -int(hpos[col]),
                                 human_aa=hrow[col], residue=aa, carriers=int(x.sum()), others=int((1 - x).sum()),
                                 families_with_both=poly, effect_log10_years=b[-1], t=t, p=2 * stats.t.sf(abs(t), df)))
    res = pd.DataFrame(rows); res["q"] = multipletests(res.p, method="fdr_bh")[1]; res = res.sort_values("p")
    res.to_csv(out, sep="\t", index=False, float_format="%.4g")
    print(f"tests: {len(res)} over {res.groupby(['gene','column']).ngroups} columns; q<0.05: {(res.q<0.05).sum()}; q<0.1: {(res.q<0.1).sum()}; smallest p {res.p.min():.2g}")
    print("p-value histogram (10 bins):", np.histogram(res.p, bins=10, range=(0, 1))[0].tolist())
    with pd.option_context("display.width", 250):
        print(res.head(25).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
