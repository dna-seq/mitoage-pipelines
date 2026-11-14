"""Inputs for Coevol: mammals spread over the tree (at most `per_family` per family, those with a body mass),
the 12 heavy-strand genes concatenated as codons, the IQ-TREE topology pruned and rooted, lifespan and mass.

Writes mammals.ali, mammals.tree and mammals.lht.
"""
import io
from pathlib import Path

import pandas as pd
from Bio import Phylo
from Bio.Data import CodonTable

from ..genes import MT_HEAVY_STRAND
from ..seqio import read_fasta


def main(data: Path, composition: Path, tree_file: Path, outgroup_file: Path | None, out: Path, per_family: int,
         quality_file: Path | None = None) -> None:
    """`outgroup_file` None: `tree_file` is the dated species tree (rooted). `quality_file`: keep only its species."""
    out.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(composition, sep="\t"); sp = sp[sp["class"] == "Mammalia"].dropna(subset=["lm"])
    tree = Phylo.read(tree_file, "newick")
    if quality_file is not None: sp = sp[sp.slug.isin(set(quality_file.read_text().split()))]
    if outgroup_file is None: sp = sp[sp.slug.isin({t.name for t in tree.get_terminals()})]
    sel = sp.sort_values("slug").groupby("family").head(per_family)
    names = {s: s.replace("-", "_") for s in sel.slug}
    if outgroup_file is not None:
        outg = outgroup_file.read_text().split()
        tree.root_with_outgroup(*[t for t in tree.get_terminals() if t.name in outg][:1])
    for t in tree.get_terminals():
        if t.name not in names: tree.prune(t)
    for cl in tree.find_clades(): cl.confidence = None; cl.branch_length = None
    for t in tree.get_terminals(): t.name = names[t.name]
    multi = [len(c.clades) for c in tree.find_clades() if len(c.clades) > 2]
    if multi: raise SystemExit(f"Coevol needs a binary tree; {len(multi)} nodes have more than two children")
    buf = io.StringIO(); Phylo.write(tree, buf, "newick", plain=True); (out / "mammals.tree").write_text(buf.getvalue())
    stops = set(CodonTable.unambiguous_dna_by_id[2].stop_codons); cat = {s: [] for s in names}
    for g in MT_HEAVY_STRAND:
        aln = read_fasta(data / "alignments" / f"{g}.codon.fasta.gz")
        rows = {s: aln[s].upper() for s in names}; L = len(next(iter(rows.values())))
        keep = [i for i in range(0, L, 3) if sum(r[i:i+3] != "---" for r in rows.values()) >= 0.9 * len(rows)]
        for s, r in rows.items():
            cod = [r[i:i+3] for i in keep]
            cat[s].append("".join("???" if (c in stops or "-" in c or any(x not in "ACGT" for x in c)) else c for c in cod))
    n = len("".join(next(iter(cat.values()))))
    (out / "mammals.ali").write_text(f"{len(cat)} {n}\n" + "".join(f"{names[s]}  {''.join(v)}\n" for s, v in cat.items()))
    (out / "mammals.lht").write_text(f"#TRAITS\n{len(sel)} 2 Longevity Mass\n" + "".join(f"{names[r.slug]}\t{r.max_longevity_years}\t{r.adult_weight_g}\n" for r in sel.itertuples()))
    print(len(sel), "species in", sel.family.nunique(), "families;", n // 3, "codons")
