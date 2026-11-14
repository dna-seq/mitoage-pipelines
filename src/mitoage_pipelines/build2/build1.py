"""MitoAge build 1 as published in the build1/ folder of the Hugging Face dataset: the 2015 CSV export,
one parquet file per CSV, values unchanged (empty cells are null). This module reads it back into what the
build-2 rules start from: one record per species, base-composition rows per region, codon rows per region.

Two facts about build 1 are rules of build 2 (site repository: pipeline/build.ts, pipeline/steps/regions.ts):
- The export has 923 species, the published release 922: Loxioides bailleui (HAGRID 1598, NC_025626.1,
  9,243 bp, not a complete mitogenome) is left out.
- The base columns are complement-labelled: column "G" holds the number of C on the RefSeq strand, "C"
  holds G, "A" holds T and "T" holds A. `fixed_counts` turns a row into true RefSeq-strand counts.
"""
import re
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from .genome import AGGREGATE_REGIONS, PROTEIN_GENES

EXCLUDED = {"Loxioides bailleui": "RefSeq NC_025626.1 is titled 'complete genome' but is only 9,243 bp "
                                  "(about 7 kb short of a complete avian mtDNA)"}
REGIONS = list(AGGREGATE_REGIONS) + list(PROTEIN_GENES)
OTHER = "Other (W\\M\\Y\\R\\S...)"
CODON = re.compile(r"^[ACGU]{3}$")


@dataclass
class Species:
    hagrid: int
    refseq: str
    kingdom: str
    phylum: str
    cls: str
    order: str
    family: str
    genus: str
    species: str
    common_name: str
    max_longevity: float


@dataclass
class Composition:
    """One base-composition row as labelled in the CSV (labels wrong, see the module docstring)."""
    label_g: int
    label_c: int
    label_a: int
    label_t: int
    other: int
    length: int
    gc_percent: float
    non_coding: int | None
    non_coding_no_dloop: int | None


@dataclass
class Build1:
    species: dict[int, Species]                       # every HAGRID of the export, 923
    kept: list[Species]                               # the 922 of the release, by HAGRID
    composition: dict[str, dict[int, list[Composition]]] = field(default_factory=dict)   # region -> hagrid -> rows
    codons: dict[str, dict[int, dict[str, int]]] = field(default_factory=dict)            # region -> hagrid -> counts


def fixed_counts(r: Composition) -> dict:
    """True RefSeq-strand counts of a build 1 row: a = column T, c = column G, g = column C, t = column A."""
    return {"a": r.label_t, "c": r.label_g, "g": r.label_c, "t": r.label_a, "other": r.other}


def _s(v) -> str:
    return "" if v is None else str(v)


def load(folder: Path) -> Build1:
    files = sorted(folder.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no build 1 parquet files in {folder}: run `mitoage download hf` first")
    species: dict[int, Species] = {}
    out = Build1(species=species, kept=[])
    for path in files:          # the CSV file names, in the order the site's loader reads them
        m = re.match(r"^(.+)_(base_composition|codon_usage)$", path.stem)
        if not m or m.group(1) not in REGIONS:
            raise ValueError(f"unexpected build 1 file {path.name}")
        region, kind = m.group(1), m.group(2)
        df = pl.read_parquet(path)
        for r in df.iter_rows(named=True):
            h = int(r["HAGRID"])
            if h not in species:
                species[h] = Species(h, r["Ref"], _s(r["Kingdom"]), _s(r["Phylum"]), _s(r["Class"]), _s(r["Order"]),
                                     _s(r["Family"]), _s(r["Genus"]), _s(r["Species"]), _s(r["Common name"]),
                                     float(r["Maximum longevity (yrs)"]))
            if kind == "base_composition":
                out.composition.setdefault(region, {}).setdefault(h, []).append(Composition(
                    int(r["G"]), int(r["C"]), int(r["A"]), int(r["T"]), int(r[OTHER]), int(r["mtDNA size (bp)"]),
                    float(r["GC%"]), r.get("length non-coding (bp)"), r.get("length non-coding without D-Loop (bp)")))
            else:
                by = out.codons.setdefault(region, {})
                if h in by:
                    raise ValueError(f"HAGRID {h} twice in {path.name}")
                counts = {c: int(r[c]) for c in df.columns if CODON.match(c)}
                if len(counts) != 64:
                    raise ValueError(f"{path.name}: expected 64 codon columns")
                by[h] = counts
    out.kept = [s for _, s in sorted(species.items()) if s.species not in EXCLUDED]
    return out
