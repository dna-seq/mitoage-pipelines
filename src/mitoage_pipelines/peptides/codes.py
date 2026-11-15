"""Genetic codes as lookup arrays over codon indices (16 a + 4 b + c, with A C G T = 0 1 2 3).

A reading of a frame is a code and a set of start codons:
  standard     NCBI table 1; start ATG only (cytosolic initiation)
  own          the species' own table (species.parquet transl_table) with its NCBI start codons; when the peptide's
               published reading is mito_agr_r, AGA and AGG are read as arginine (they are stops in table 2)
The published readings of the catalogue map onto these: "standard" is `standard`; "mito" and "mito_agr_r" are
`own` for a species of table 2 (the human's).
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from Bio.Data import CodonTable

NT = "ACGT"
CODONS = [a + b + c for a in NT for b in NT for c in NT]
AGR = ("AGA", "AGG")
INDEX = np.full(256, 4, dtype=np.int8)        # byte -> 0..3, anything else (N, gaps, IUPAC) -> 4
for _i, _c in enumerate(NT):
    INDEX[ord(_c)] = _i


@dataclass(frozen=True)
class Table:
    """One reading: per codon index, the residue, whether it stops, whether it starts."""
    aa: tuple[str, ...]
    stop: np.ndarray
    start: np.ndarray

    def is_stop(self, codon: str) -> bool:
        return codon in CODONS and bool(self.stop[CODONS.index(codon)])

    def is_start(self, codon: str) -> bool:
        return codon in CODONS and bool(self.start[CODONS.index(codon)])


@lru_cache(maxsize=None)
def tables(table_id: int, reading: str) -> Table:
    """The reading `reading` (standard, mito, mito_agr_r, own, own_agr_r) for NCBI table `table_id`. standard ignores
    the table. mito/own are the table as NCBI has it; *_agr_r reads AGA and AGG as arginine."""
    if reading == "standard":
        t = CodonTable.unambiguous_dna_by_id[1]
        starts = {"ATG"}
    else:
        t = CodonTable.unambiguous_dna_by_id[table_id]
        starts = set(t.start_codons)
    fwd = dict(t.forward_table)
    stops = set(t.stop_codons)
    if reading.endswith("agr_r"):
        for c in AGR:
            if c in stops:
                stops.discard(c); fwd[c] = "R"
    aa = tuple("*" if c in stops else fwd[c] for c in CODONS)
    return Table(aa, np.array([c in stops for c in CODONS]), np.array([c in starts for c in CODONS]))


def species_reading(reading: str, published: str) -> str:
    """The reading name `tables` takes for one of this analysis' two readings (own, standard) of a peptide whose
    published reading is `published`."""
    if reading == "standard":
        return "standard"
    return "own_agr_r" if published == "mito_agr_r" else "own"


def translate(nt: str, table: Table) -> str:
    """Codon by codon; a codon with a base other than ACGT is X."""
    out = []
    for i in range(0, len(nt) - 2, 3):
        c = nt[i:i + 3]
        out.append(table.aa[CODONS.index(c)] if c in CODONS else "X")
    return "".join(out)


def to_index(seq: str) -> np.ndarray:
    """Bases as 0..3, anything else 4."""
    return INDEX[np.frombuffer(seq.encode(), dtype=np.uint8)]


def codon_index(nt: np.ndarray) -> np.ndarray:
    """Index of the codon starting at each position of `nt` (0..3 per base, 4 unknown), -1 when a base is unknown;
    the last two positions, which start no complete codon, are -2."""
    n = len(nt)
    out = np.full(n, -2, dtype=np.int16)
    if n >= 3:
        a, b, c = nt[:-2].astype(np.int16), nt[1:-1].astype(np.int16), nt[2:].astype(np.int16)
        idx = 16 * a + 4 * b + c
        idx[(a > 3) | (b > 3) | (c > 3)] = -1
        out[:n - 2] = idx
    return out
