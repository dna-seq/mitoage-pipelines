"""NCBI Taxonomy (new_taxdump): ranks and parents, merged taxids, ranked lineages and names, as build 2
reads them (site repository: pipeline/sources/ncbi-taxonomy.ts).

- The current taxid of a merged one follows merged.dmp.
- The species of a taxon is its nearest ancestor (or itself) of rank species, inside Metazoa only.
- Names match on a key: parenthesised parts and quotes removed, single spaces, lower case. Scientific
  names and synonyms (synonym, equivalent name, genbank synonym) are separate indexes, each name key
  mapping to the sorted species taxids of the Metazoa taxa at or below species rank that carry it.
"""
import re
from pathlib import Path

import numpy as np
import polars as pl

METAZOA = 33208
SYN_CLASSES = {"synonym", "equivalent name", "genbank synonym"}


def name_key(s: str) -> str:
    s = re.sub(r"\([^)]*\)", " ", s)
    s = re.sub("[\u2018\u2019'\"]", " ", s)
    return re.sub(r"\s+", " ", s).strip().lower()


def _dmp(path: Path, cols: list[int], names: list[str]) -> pl.DataFrame:
    """A .dmp file (fields separated by tab-bar-tab): the given field positions."""
    if not path.is_file():
        raise FileNotFoundError(f"no {path}: run `mitoage download taxonomy` first")
    df = pl.read_csv(path, separator="\t", has_header=False, quote_char=None, infer_schema=False,
                     truncate_ragged_lines=True, columns=[2 * c for c in cols])
    return df.rename(dict(zip(df.columns, names)))


class Taxonomy:
    def __init__(self, folder: Path):
        self.folder = folder
        nodes = _dmp(folder / "nodes.dmp", [0, 1, 2], ["taxid", "parent", "rank"])
        tax = nodes["taxid"].cast(pl.Int64).to_numpy(); par = nodes["parent"].cast(pl.Int64).to_numpy()
        size = int(tax.max()) + 1
        self.parent = np.zeros(size, dtype=np.int64); self.parent[tax] = par
        rank = nodes["rank"].to_numpy()
        self.is_species = np.zeros(size, dtype=bool); self.is_species[tax[rank == "species"]] = True
        self.rank_of = dict(zip(tax[(rank != "species") & (rank != "no rank")].tolist(), rank[(rank != "species") & (rank != "no rank")].tolist()))
        merged = _dmp(folder / "merged.dmp", [0, 1], ["old", "new"])
        self.merged = dict(zip(merged["old"].cast(pl.Int64).to_list(), merged["new"].cast(pl.Int64).to_list()))
        self._meta: dict[int, bool] = {}
        self._sp: dict[int, int] = {}
        self._lineage: dict[int, dict] | None = None
        self.sci: dict[str, list[int]] = {}
        self.syn: dict[str, list[int]] = {}
        self.names: dict[int, dict] = {}
        self.label = folder.name

    def exists(self, t: int) -> bool:
        return 0 < t < len(self.parent) and self.parent[t] != 0

    def current(self, t: int) -> int:
        return self.merged.get(t, t)

    def is_metazoa(self, t: int) -> bool:
        if not self.exists(t): return False
        path, x = [], t
        while True:
            if x in self._meta: res = self._meta[x]; break
            path.append(x)
            if x == METAZOA: res = True; break
            p = int(self.parent[x])
            if p == x or p == 0: res = False; break
            x = p
        for y in path: self._meta[y] = res
        return res

    def species_of(self, t0: int) -> int | None:
        """Nearest species-rank ancestor (or self) of a Metazoa taxon; None above species or outside Metazoa."""
        t = self.current(t0)
        if not self.is_metazoa(t): return None
        path, x = [], t
        while True:
            if x in self._sp: found = self._sp[x]; break
            path.append(x)
            if self.is_species[x]: found = x; break
            p = int(self.parent[x])
            if x == METAZOA or p == x or p == 0: found = 0; break
            x = p
        for y in path: self._sp[y] = found
        return found or None

    def lineage(self, t: int) -> dict | None:
        """name, species, genus, family, order, class, phylum, kingdom of a taxid (rankedlineage.dmp)."""
        if self._lineage is None:
            cols = ["taxid", "name", "species", "genus", "family", "order", "class", "phylum", "kingdom"]
            df = _dmp(self.folder / "rankedlineage.dmp", list(range(9)), cols).with_columns(pl.col("taxid").cast(pl.Int64))
            df = df.with_columns([pl.col(c).fill_null("") for c in cols[1:]])
            self._lineage = {r[0]: dict(zip(cols[1:], r[1:])) for r in df.iter_rows()}
        return self._lineage.get(t)

    def load_names(self, keys: set[str], taxids: set[int]) -> None:
        """The matching indexes (restricted to `keys`) and the full name lists of `taxids`."""
        df = _dmp(self.folder / "names.dmp", [0, 1, 3], ["taxid", "name", "cls"]).with_columns(pl.col("taxid").cast(pl.Int64))
        df = df.filter(pl.col("cls").is_in(["scientific name", "genbank common name", *SYN_CLASSES]))
        for t, name, cls in df.filter(pl.col("taxid").is_in(list(taxids))).iter_rows():
            n = self.names.setdefault(t, {"scientific": None, "synonyms": [], "genbank_common": None})
            if cls == "scientific name": n["scientific"] = name
            elif cls in SYN_CLASSES: n["synonyms"].append(name)
            elif cls == "genbank common name" and n["genbank_common"] is None: n["genbank_common"] = name
        key = (pl.col("name").str.replace_all(r"\([^)]*\)", " ").str.replace_all("[\u2018\u2019'\"]", " ")
               .str.replace_all(r"\s+", " ").str.strip_chars().str.to_lowercase())
        df = df.filter(pl.col("cls") != "genbank common name").with_columns(key.alias("key")).filter(pl.col("key").is_in(list(keys)))
        sci: dict[str, set[int]] = {}; syn: dict[str, set[int]] = {}
        for t, k, cls in df.select("taxid", "key", "cls").iter_rows():
            sp = self.species_of(t)
            if sp is None: continue
            (sci if cls == "scientific name" else syn).setdefault(k, set()).add(sp)
        self.sci = {k: sorted(v) for k, v in sci.items()}
        self.syn = {k: sorted(v) for k, v in syn.items()}
