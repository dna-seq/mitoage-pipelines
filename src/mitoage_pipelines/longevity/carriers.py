"""Which mammalian families carry the minor residue at the two columns the phylogenetic test puts first
(ND5 column 316, residue H; CYTB column 83, residue M), and how long those species live. Prints only."""
from pathlib import Path

import pandas as pd

from ..seqio import read_fasta

# gene, alignment column (1-based), the residue whose carriers are listed
SITES = (("ND5", 316, "H"), ("CYTB", 83, "M"))


def main(data: Path, composition: Path) -> None:
    sp = pd.read_csv(composition, sep="\t"); sp = sp[sp["class"] == "Mammalia"].dropna(subset=["lm"])
    for gene, col, residue in SITES:
        rows = read_fasta(data / "alignments" / f"{gene}.aa.fasta.gz")
        sp["res"] = sp.slug.map(lambda s: rows[s][col - 1])
        print(gene, col, sp.res.value_counts().to_dict())
        minor = sp[sp.res.isin([residue])]
        print(minor.groupby(["order", "family"]).agg(n=("slug", "size"), median_years=("max_longevity_years", "median")).to_string())
        print("  family totals:", sp[sp.family.isin(minor.family)].groupby("family").size().to_dict())
