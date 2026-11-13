"""How many of the build's mammals each NCBI ortholog set covers. Prints only."""
import re
from pathlib import Path

import polars as pl

from ..seqio import read_records


def gene_dirs(ncbi: Path) -> list[str]:
    genes = sorted(p.name for p in ncbi.iterdir() if p.is_dir()) if ncbi.is_dir() else []
    if not genes:
        raise FileNotFoundError(f"no ortholog sets under {ncbi}: run `mitoage download orthologs` first")
    return genes


def main(data: Path, ncbi: Path) -> None:
    sp = pl.read_parquet(data / "species.parquet").filter(pl.col("class") == "Mammalia")
    ours = {" ".join(n.split()[:2]) for n in sp["scientific_name"]}
    withlife = {" ".join(n.split()[:2]) for n in sp.filter(pl.col("max_longevity_years") > 0)["scientific_name"]}
    genes = gene_dirs(ncbi); sets = {}
    for g in genes:
        s = set()
        for _, description, _ in read_records(ncbi / g / "ncbi_dataset" / "data" / "protein.faa"):
            m = re.search(r"\[organism=([^\]]+)\]", description); s.add(" ".join(m.group(1).split()[:2]))
        sets[g] = s; print(f"{g}\tspecies at NCBI {len(s)}\tin our mammals {len(s & ours)}")
    allg = set.intersection(*sets.values()); n8 = {x for x in ours if sum(x in s for s in sets.values()) >= 9}
    print(f"our mammals: {len(ours)} (with lifespan {len(withlife)}); with all {len(genes)} subunits: {len(allg & ours)}; with 9 or more: {len(n8)}")
    print("union at NCBI:", len(set.union(*sets.values())), "; our mammals with at least one:", len(set.union(*sets.values()) & ours))
