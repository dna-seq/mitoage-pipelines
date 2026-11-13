"""Do nuclear subunits of the respiratory complexes change rate on the same branches as the
mitochondrial genome, more than nuclear proteins that do not touch mitochondrial ones?

Branch lengths per protein on one topology (from `branch-lengths`). Per branch:
  nuclear background = mean length over the cytosolic control proteins (each scaled to mean 1)
  mt rate            = mean length over the 13 mt proteins (scaled likewise) / nuclear background
  a nuclear protein's relative rate = its scaled length / nuclear background
Test A: Spearman correlation over branches of a protein's relative rate with the mt rate, per group.
Test B: the same against one complex's mt subunits only, the mt genes scaled by the mt mean.
Leave-one-out: a control protein's own length is left out of its background.

Writes erc.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats

from .. import genes as gene_lists


def branches(path):
    """bipartition (its smaller side, as a set of tip names) -> branch length"""
    t = Phylo.read(path, "newick"); alltips = frozenset(x.name for x in t.get_terminals()); out = {}
    for cl in t.find_clades():
        if cl is t.root or cl.branch_length is None: continue
        s = frozenset(x.name for x in cl.get_terminals()); key = min(s, alltips - s, key=lambda z: (len(z), sorted(z)))
        out[key] = out.get(key, 0.0) + cl.branch_length
    return out


def main(genes_file: Path, bl_dir: Path, out: Path) -> None:
    """`out` is the table to write."""
    out.parent.mkdir(parents=True, exist_ok=True)
    GROUPS = {k: list(v) for k, v in gene_lists.GROUPS.items()}
    genes = genes_file.read_text().split()
    known = {g for v in GROUPS.values() for g in v}
    GROUPS["cytosolic control"] = [g for g in genes if not g.startswith("MT-") and g not in known]
    MTC = gene_lists.MT_OF_COMPLEX
    bl = {g: branches(bl_dir / f"{g}.treefile") for g in genes if (bl_dir / f"{g}.treefile").exists()}
    keys = sorted(set.intersection(*[set(b) for b in bl.values()]), key=lambda z: (len(z), sorted(z)))
    M = pd.DataFrame({g: [b[k] for k in keys] for g, b in bl.items()}); M = M / M.mean()
    term = np.array([len(k) == 1 for k in keys])
    ctrl = [g for g in GROUPS["cytosolic control"] if g in M]; mt = [g for g in M if g.startswith("MT-")]
    bg = M[ctrl].mean(axis=1); mtmean = M[mt].mean(axis=1)
    ok = (bg > np.quantile(bg, 0.1)) & (mtmean > 0)             # the shortest tenth of branches carry no rate information
    print(f"proteins {M.shape[1]} (mt {len(mt)}, cytosolic controls {len(ctrl)}); branches {len(keys)}, used {int(ok.sum())} ({int((ok & term).sum())} terminal)")
    print(f"mt rate against nuclear background over branches: Spearman {stats.spearmanr(mtmean[ok], bg[ok])[0]:.2f}")
    rows = []
    for grp, gs in GROUPS.items():
        for g in gs:
            if g not in M: continue
            b = M[[c for c in ctrl if c != g]].mean(axis=1) if g in ctrl else bg
            rel = M[g] / b; mtrate = mtmean / b
            r = {"group": grp, "gene": g, "A_mt_rate": stats.spearmanr(rel[ok], mtrate[ok])[0]}
            for cx, ms in MTC.items():
                r[f"B_complex_{cx}"] = stats.spearmanr(rel[ok], (M[ms].mean(axis=1) / mtmean)[ok])[0]
            rows.append(r)
    res = pd.DataFrame(rows); res.to_csv(out, sep="\t", index=False, float_format="%.3f")
    c = res[res.group == "cytosolic control"]
    print("\nper group: median correlation; Mann-Whitney p against the cytosolic controls")
    for grp in GROUPS:
        d = res[res.group == grp]
        if d.empty: continue
        line = f"{grp:28s} n={len(d):2d}"
        for col in ["A_mt_rate", "B_complex_IV", "B_complex_III"]:
            p = stats.mannwhitneyu(d[col], c[col], alternative="two-sided")[1] if grp != "cytosolic control" and len(d) > 1 else float("nan")
            line += f"  {col} {d[col].median():+.2f} (p {p:.2g})"
        print(line)
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(res.sort_values("A_mt_rate", ascending=False).to_string(index=False, float_format=lambda v: f"{v:+.2f}"))
