"""Preparation for HyPhy RELAX: long-lived against short-lived mammals.

Lifespan residual = log10 lifespan minus what body mass predicts (plain regression within mammals).
Foreground and Reference are balanced within orders (see below), so that the contrast is not one
order against another. The tree is the IQ-TREE one, pruned to these species, tips labelled
{Foreground} / {Reference}.

Writes species.tsv, tree.nwk and one codon alignment per gene (<GENE>.fna).
"""
import io
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Data import CodonTable

from ..genes import MT_GENES
from ..seqio import read_fasta


def main(data: Path, composition: Path, tree_file: Path, out: Path, n: int, per_family: int,
         exclude: Path | None = None) -> None:
    """`n`: species per order at each end at most; `per_family`: of one family at each end at most;
    `exclude`: the species.tsv of an earlier draw, whose species are left out (a second, disjoint draw)."""
    out.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(composition, sep="\t"); sp = sp[sp["class"] == "Mammalia"].dropna(subset=["lm"])
    b = np.polyfit(sp.lm, sp.ll, 1); sp["resid"] = sp.ll - np.polyval(b, sp.lm)
    if exclude is not None:      # a second, disjoint draw: leave out the species of an earlier one
        sp = sp[~sp.slug.isin(pd.read_csv(exclude, sep="\t").slug)]
    # balanced within orders: in every order with 10 or more species, the same number from each end,
    # a fifth of the order up to n, at most per_family of a family at each end
    fgs, refs = [], []
    for o, d in sp.groupby("order"):
        if len(d) < 10: continue
        k = min(n, len(d) // 5)
        hi = d.sort_values("resid", ascending=False).groupby("family").head(per_family).head(k)
        lo = d[~d.slug.isin(hi.slug)].sort_values("resid").groupby("family").head(per_family).head(k)
        k = min(len(hi), len(lo)); fgs.append(hi.head(k)); refs.append(lo.head(k))
    fg, ref = pd.concat(fgs), pd.concat(refs)
    sel = pd.concat([fg.assign(set="Foreground"), ref.assign(set="Reference")])
    sel[["slug", "order", "family", "max_longevity_years", "adult_weight_g", "resid", "set"]].to_csv(out / "species.tsv", sep="\t", index=False, float_format="%.4g")
    print(sel.groupby("set").agg(n=("slug", "size"), families=("family", "nunique"), median_years=("max_longevity_years", "median"), median_g=("adult_weight_g", "median"), resid=("resid", "median")).to_string())
    print("orders:", sel.groupby(["set", "order"]).size().unstack(0).fillna(0).astype(int).to_string())
    label = dict(zip(sel.slug, sel.set))
    tree = Phylo.read(tree_file, "newick")
    for t in tree.get_terminals():
        if t.name not in label: tree.prune(t)
    for cl in tree.find_clades(): cl.confidence = None
    for t in tree.get_terminals(): t.name = t.name.replace("-", "_") + "{" + label[t.name] + "}"
    buf = io.StringIO(); Phylo.write(tree, buf, "newick"); nwk = buf.getvalue().replace("'", "")
    (out / "tree.nwk").write_text(nwk)
    stops = set(CodonTable.unambiguous_dna_by_id[2].stop_codons)
    for g in MT_GENES:
        aln = read_fasta(data / "alignments" / f"{g}.codon.fasta.gz")
        rows = {s: aln[s].upper() for s in label}
        L = len(next(iter(rows.values()))); assert L % 3 == 0
        keep = [i for i in range(0, L, 3) if sum(r[i:i+3] != "---" for r in rows.values()) >= 0.5 * len(rows)]
        masked = 0
        with open(out / f"{g}.fna", "w") as f:
            for s, r in rows.items():
                cod = [r[i:i+3] for i in keep]
                for j, c in enumerate(cod):
                    if c in stops or ("-" in c and c != "---") or any(x not in "ACGT-" for x in c): cod[j] = "---"; masked += c != "---"
                f.write(f">{s.replace('-', '_')}\n{''.join(cod)}\n")
        print(g, "codons kept", len(keep), "codons masked (stops, partial, ambiguous)", masked)
