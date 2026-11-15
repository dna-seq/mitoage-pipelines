"""Compare the result tables of two runs of slurm/manuscript-rerun.sbatch, table by table and class by class.

    uv run python slurm/compare-runs.py <root A> <root B>

Each argument is a data root (MITOAGE_PIPELINES_DATA); the tables are read from <root>/output/longevity. For each
table in either run it prints whether the two files are byte-identical, and otherwise, for each value of the
`class` column, how many rows each run has and how many rows of A are found unchanged in B (the whole row, as
text). A class whose rows are all shared is "identical"; a table without a class column is compared as a whole.
A table present in one run only is reported as such, never skipped.
"""
import sys
from collections import defaultdict
from pathlib import Path


def rows_by_class(path: Path) -> tuple[str, dict[str, list[str]]]:
    lines = path.read_text().splitlines()
    header, body = lines[0], lines[1:]
    cols = header.split("\t")
    by: dict[str, list[str]] = defaultdict(list)
    if "class" in cols:
        i = cols.index("class")
        for line in body:
            by[line.split("\t")[i]].append(line)
    else:
        by["(whole table)"] = body
    return header, by


def main(a: Path, b: Path) -> None:
    da, db = a / "output" / "longevity", b / "output" / "longevity"
    names = sorted({p.name for p in da.glob("*.tsv")} | {p.name for p in db.glob("*.tsv")})
    for name in names:
        pa, pb = da / name, db / name
        if not pa.is_file() or not pb.is_file():
            print(f"{name}: only in {'A' if pa.is_file() else 'B'}")
            continue
        if pa.read_bytes() == pb.read_bytes():
            print(f"{name}: identical bytes")
            continue
        ha, ra = rows_by_class(pa)
        hb, rb = rows_by_class(pb)
        if ha != hb:
            print(f"{name}: the headers differ")
            continue
        print(f"{name}: differs")
        for cls in sorted(set(ra) | set(rb)):
            xa, xb = ra.get(cls, []), rb.get(cls, [])
            shared = len(set(xa) & set(xb))
            state = "identical" if sorted(xa) == sorted(xb) else "differs"
            print(f"  {cls}: {state}; rows A {len(xa)}, B {len(xb)}, shared {shared}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(Path(sys.argv[1]), Path(sys.argv[2]))
