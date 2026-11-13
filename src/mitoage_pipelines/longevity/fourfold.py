"""Base composition at fourfold-degenerate third positions, per species, over the 12 heavy-strand genes.

Writes fourfold.tsv: slug, the shares of A, C and T among the fourfold sites, and their number.
"""
from pathlib import Path

from ..genes import MT_HEAVY_STRAND
from ..seqio import read_fasta

# first two bases of the codon families that are fourfold degenerate in every animal mitochondrial code
FOUR = {"CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"}


def count_fourfold(seqs: dict[str, str], cnt: dict[str, dict[str, int]] | None = None) -> dict[str, dict[str, int]]:
    """Add to `cnt` (slug -> base -> count) the third bases of the fourfold codons of a codon alignment.

    Codons are read in frame from the first column; a third position that is not A, C, G or T is not counted.
    """
    cnt = {} if cnt is None else cnt
    for slug, s in seqs.items():
        s = s.upper()
        c = cnt.setdefault(slug, dict.fromkeys("ACGT", 0))
        for i in range(0, len(s) - 2, 3):
            if s[i:i+2] in FOUR and s[i+2] in c: c[s[i+2]] += 1
    return cnt


def fourfold_counts(data: Path) -> dict[str, dict[str, int]]:
    cnt: dict[str, dict[str, int]] = {}
    for g in MT_HEAVY_STRAND:
        count_fourfold(read_fasta(data / "alignments" / f"{g}.codon.fasta.gz"), cnt)
    return cnt


def main(data: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    cnt = fourfold_counts(data)
    with open(out / "fourfold.tsv", "w") as f:
        f.write("slug\tA4\tC4\tT4\tfourfold_sites\n")
        for s, c in sorted(cnt.items()):
            n = max(1, sum(c.values())); f.write(f"{s}\t{c['A']/n:.5f}\t{c['C']/n:.5f}\t{c['T']/n:.5f}\t{n}\n")
    print(f"fourfold.tsv: {len(cnt)} species")
