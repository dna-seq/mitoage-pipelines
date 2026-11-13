"""The one FASTA reader, on polars-bio.

`polars_bio.read_fasta` returns the columns `name` (the header up to the first whitespace),
`description` (the rest of the header, null when there is none) and `sequence` (the lines joined,
case and gap characters kept). It reads gzip by the `.gz` extension, and plain FASTA under any
extension tried here (.fasta, .fa, .faa, .fna, .afa), so no extension needs a second parser.
"""
import os
from pathlib import Path

os.environ.setdefault("TQDM_DISABLE", "1")     # polars-bio draws a progress bar per file on stderr

import polars_bio as pb  # noqa: E402


def read_records(path: Path) -> list[tuple[str, str, str]]:
    """Every record in file order: (name, description, sequence); the description is "" when absent."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"no FASTA file at {path}")
    df = pb.read_fasta(str(path))
    return [(n, d or "", s or "") for n, d, s in zip(df["name"], df["description"], df["sequence"])]


def read_fasta(path: Path) -> dict[str, str]:
    """Sequences keyed by the first whitespace-delimited token of the header, in file order.

    Sequences are returned as written (case and gaps kept): each caller upper-cases or strips gaps
    itself. A name that occurs twice is refused: a dictionary would keep one of the two silently.
    """
    out: dict[str, str] = {}
    for name, _, seq in read_records(path):
        if name in out:
            raise ValueError(f"{path}: duplicate sequence name {name}")
        out[name] = seq
    return out


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    Path(path).write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
