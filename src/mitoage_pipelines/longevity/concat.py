"""Concatenate the 13 protein alignments for every species of one class that has all of them,
plus six species from each of two other classes as outgroup (of Aves, Reptilia, Mammalia: the two
that are not the class). Columns under half occupied are dropped, per gene.

Writes <Class>.faa and <Class>.outgroup.txt.
"""
from pathlib import Path

import pandas as pd

from ..genes import MT_GENES
from ..seqio import read_fasta


def main(data: Path, composition: Path, out: Path, cls: str = "Mammalia") -> None:
    out.mkdir(parents=True, exist_ok=True)
    og = [c for c in ("Aves", "Reptilia", "Mammalia") if c != cls][:2]
    sp = pd.read_csv(composition, sep="\t")          # species with all 13 proteins and a lifespan
    mam = sorted(sp[sp["class"] == cls].slug)
    outg = sorted(sp[sp["class"] == og[0]].slug)[:6] + sorted(sp[sp["class"] == og[1]].slug)[:6]
    seqs = {s: [] for s in mam + outg}; parts = []; pos = 1
    for g in MT_GENES:
        rows = {k: v for k, v in read_fasta(data / "alignments" / f"{g}.aa.fasta.gz").items() if k in seqs}
        keep = [i for i in range(len(next(iter(rows.values())))) if sum(r[i] != "-" for r in rows.values()) >= 0.5 * len(rows)]
        for s in seqs: seqs[s].append("".join(rows[s][i] for i in keep))
        parts.append(f"AA, {g} = {pos}-{pos + len(keep) - 1}"); pos += len(keep)
    (out / f"{cls}.faa").write_text("".join(f">{s}\n{''.join(v)}\n" for s, v in seqs.items()))
    (out / f"{cls}.outgroup.txt").write_text("\n".join(outg) + "\n")
    print(len(mam), cls + ",", len(outg), "outgroup,", pos - 1, "columns")
