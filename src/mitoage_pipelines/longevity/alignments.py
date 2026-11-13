"""Which protein alignment a per-column test reads: the build's, or another tool's of the same sequences."""
from pathlib import Path

from ..seqio import read_fasta


def protein_alignment(data: Path, gene: str, alt_dir: Path | None = None, alt_tool: str | None = None) -> dict[str, str]:
    """The build's `alignments/<GENE>.aa.fasta.gz`, or `<alt_dir>/<GENE>.<alt_tool>.afa` when both are given."""
    if (alt_dir is None) != (alt_tool is None):
        raise ValueError("an alternative alignment needs both its folder and the tool's name")
    if alt_dir is not None:
        return read_fasta(alt_dir / f"{gene}.{alt_tool}.afa")
    return read_fasta(data / "alignments" / f"{gene}.aa.fasta.gz")
