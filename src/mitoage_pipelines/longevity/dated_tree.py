"""The build's dated species tree as newick with species slugs, and the species whose lifespan record is of
high or acceptable quality. Both were fixed as the analysis' ground before any result: the tree is never one
made from the genes under test."""
import json
import sys
from pathlib import Path

import polars as pl

from . import species_sets

QUALITY_OK = ("high", "acceptable")


def to_newick(node: dict, name_of: dict) -> str:
    """One node of phylogeny.json (`length`, and either `children` or the species id) as newick."""
    children = node.get("children")
    length = node.get("length", 0)
    if children:
        return "(" + ",".join(to_newick(c, name_of) for c in children) + f"):{length}"
    key = next(k for k in node if k not in ("length", "children"))
    return f"{name_of[node[key]]}:{length}"


def main(data: Path, site_tree: Path, tree_out: Path, quality_out: Path) -> None:
    sp = pl.read_parquet(data / "species.parquet")
    doc = json.loads(site_tree.read_text())
    sys.setrecursionlimit(100000)
    newick = to_newick(doc["root"], dict(zip(sp["id"], sp["slug"])))
    tree_out.parent.mkdir(parents=True, exist_ok=True)
    tree_out.write_text(newick + ";\n")
    ok = sp.filter(pl.col("lifespan_quality").is_in(QUALITY_OK))
    name = species_sets.active()
    if name is not None:        # H27: the analyses that take their species from this list keep the set's only
        n0 = len(ok)
        ok = ok.filter(pl.col("slug").is_in(sorted(species_sets.slugs(data, name))))
        print(f"species set {name}: {len(ok)} of the {n0} species with a usable lifespan record")
    quality_out.write_text("\n".join(ok["slug"]) + "\n")
    print(f"dated tree {doc['version']}: {newick.count(',') + 1} species; lifespan quality high or acceptable: {len(ok)} of {len(sp)}")
