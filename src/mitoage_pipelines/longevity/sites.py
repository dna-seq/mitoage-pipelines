"""Every alignment column of the 13 proteins against lifespan, in mammals, without a tree.

Per column and residue: is carrying this residue associated with log10 maximum lifespan, comparing only
within families, with body mass and fourfold base composition in the model. A test needs the residue and
its absence in at least 15 species each, and at least 5 families that hold both.

Kept as the cautionary version: families are not independent, and the test is inflated; `pgls` is the
test with a tree.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..genes import MT_GENES
from .alignments import protein_alignment


def centre(v, f):
    m = np.zeros(f.max() + 1); n = np.bincount(f, minlength=f.max() + 1)
    np.add.at(m, f, v); return v - (m / np.maximum(n, 1))[f]


def main(data: Path, composition: Path, fourfold: Path, out: Path,
         alt_dir: Path | None = None, alt_tool: str | None = None) -> None:
    """`out` is the table to write. `alt_dir`/`alt_tool`: read <alt_dir>/<GENE>.<alt_tool>.afa instead of the build's alignments."""
    out.parent.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(composition, sep="\t")
    comp = pd.read_csv(fourfold, sep="\t")
    sp = sp[sp["class"] == "Mammalia"].dropna(subset=["lm"]).merge(comp, on="slug")
    sp = sp.groupby("family").filter(lambda x: len(x) >= 3).reset_index(drop=True)
    print(f"mammals with mass in families of 3 or more: {len(sp)} in {sp.family.nunique()} families")
    fam = sp.family.astype("category").cat.codes.to_numpy()
    Y = sp.ll.to_numpy(); COV = sp[["lm", "A4", "C4", "T4"]].to_numpy()

    rows = []
    for g in MT_GENES:
        aln = protein_alignment(data, g, alt_dir, alt_tool)
        M = np.array([list(aln[s]) for s in sp.slug])
        hrow = aln["homo-sapiens"]; hpos = np.cumsum([c != "-" for c in hrow])
        for col in range(M.shape[1]):
            c = M[:, col]; ok = c != "-"
            if ok.mean() < 0.8: continue
            for aa in np.unique(c[ok]):
                x = (c[ok] == aa).astype(float)
                if x.sum() < 15 or (1 - x).sum() < 15: continue
                f = fam[ok]; _, f = np.unique(f, return_inverse=True)
                nf = np.bincount(f); s1 = np.bincount(f, weights=x)
                poly = int(((s1 > 0) & (s1 < nf)).sum())
                if poly < 5: continue
                X = np.column_stack([centre(x, f)] + [centre(COV[ok, j], f) for j in range(COV.shape[1])])
                y = centre(Y[ok], f)
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                resid = y - X @ beta; df = len(y) - len(nf) - X.shape[1]
                cov = np.linalg.pinv(X.T @ X) * (resid @ resid / df)
                t = beta[0] / np.sqrt(cov[0, 0])
                rows.append(dict(gene=g, column=col + 1, human_pos=int(hpos[col]) if hrow[col] != "-" else -int(hpos[col]),
                                 human_aa=hrow[col], residue=aa, carriers=int(x.sum()), others=int((1 - x).sum()),
                                 families_with_both=poly, effect_log10_years=beta[0], t=t, p=2 * stats.t.sf(abs(t), df)))
    res = pd.DataFrame(rows)
    res["q"] = multipletests(res.p, method="fdr_bh")[1]
    res = res.sort_values("p")
    res.to_csv(out, sep="\t", index=False, float_format="%.4g")
    print(f"tests: {len(res)} over {res.groupby(['gene','column']).ngroups} columns; q<0.05: {(res.q<0.05).sum()}; q<0.1: {(res.q<0.1).sum()}; smallest p {res.p.min():.2g}; expected smallest by chance about {1/len(res):.2g}")
    print("p-value histogram (10 bins):", np.histogram(res.p, bins=10, range=(0, 1))[0].tolist())
    with pd.option_context("display.width", 250):
        print(res.head(30).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
        print(res.groupby("gene").agg(tests=("p", "size"), q10=("q", lambda q: (q < 0.1).sum())).to_string())
