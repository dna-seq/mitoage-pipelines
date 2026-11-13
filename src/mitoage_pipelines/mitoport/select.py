"""One protein and its coding sequence per species for each nuclear gene, from NCBI's ortholog sets.
Per species: the isoform whose length is nearest the human reference's; ties to the first listed.
A species with two gene ids in one set is counted, and the gene nearer the human length kept.

Writes <GENE>.faa and <GENE>.fna, headers `<species slug> <protein accession>`.
"""
import re
from pathlib import Path

import polars as pl

from ..seqio import read_records
from .overlap import gene_dirs

B = "TCAG"; AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}


def tr(c):
    return "".join(CODE.get(c[i:i+3].upper(), "X") for i in range(0, len(c) - 2, 3)).rstrip("*")


def fasta(p: Path):
    """(whole header, sequence) of each record, in file order."""
    for name, description, seq in read_records(p):
        yield f"{name} {description}".strip(), seq


def main(data: Path, ncbi: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    sp = pl.read_parquet(data / "species.parquet").filter(pl.col("class") == "Mammalia")
    slug = {" ".join(n.split()[:2]): s for n, s in zip(sp["scientific_name"], sp["slug"])}
    for g in gene_dirs(ncbi):
        d = ncbi / g / "ncbi_dataset" / "data"
        # the report does not link transcript to protein: match each CDS to a protein of the same gene by translation
        bygene = {}
        for h, c in fasta(d / "cds.fna"):
            bygene.setdefault(re.search(r"GeneID=(\d+)", h).group(1), []).append(c)
        proteins = list(fasta(d / "protein.faa"))
        cds = {}
        for h, s in proteins:
            gid = re.search(r"GeneID=(\d+)", h).group(1)
            for c in bygene.get(gid, []):
                if tr(c).replace("U", "*") == s or tr(c) == s: cds[h.split()[0]] = c; break
        prots = [(h.split()[0], re.search(r"\[organism=([^\]]+)\]", h).group(1), re.search(r"GeneID=(\d+)", h).group(1), s) for h, s in proteins]
        hlen = next(len(s) for a, o, gid, s in prots if o == "Homo sapiens" and a.startswith("NP_"))
        best, genes = {}, {}
        for a, o, gid, s in prots:
            k = " ".join(o.split()[:2]); genes.setdefault(k, set()).add(gid)
            if k not in slug or a not in cds: continue
            if k not in best or abs(len(s) - hlen) < abs(len(best[k][2]) - hlen): best[k] = (a, gid, s)
        multi = {k: v for k, v in genes.items() if len(v) > 1 and k in slug}
        bad = 0
        with open(out / f"{g}.faa", "w") as fa, open(out / f"{g}.fna", "w") as fn:
            for k, (a, gid, s) in sorted(best.items()):
                c = cds[a]
                if len(c) != 3 * (len(s) + 1) and len(c) != 3 * len(s): bad += 1; continue
                fa.write(f">{slug[k]} {a}\n{s}\n"); fn.write(f">{slug[k]} {a}\n{c[:3 * len(s)]}\n")
        lens = sorted(len(v[2]) for v in best.values())
        print(f"{g}\thuman {hlen} aa\tspecies {len(best)}\tlength min/median/max {lens[0]}/{lens[len(lens)//2]}/{lens[-1]}\toff by >20% {sum(abs(l-hlen) > 0.2*hlen for l in lens)}\tspecies with 2+ gene ids {len(multi)}\tCDS length mismatch dropped {bad}")
