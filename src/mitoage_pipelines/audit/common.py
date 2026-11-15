"""Pieces every audit step shares: a three-valued answer, warnings counted by reason, and tables written in a
stable order with fixed number formats, so that a rerun gives the same bytes."""
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

YES, NO, UNKNOWN = "yes", "no", "unknown"


class Warnings:
    """Warnings raised in a loop, kept by reason: one line per reason with its count and the first example."""

    def __init__(self, step: str):
        self.step = step
        self.count: Counter = Counter()
        self.example: dict[str, str] = {}

    def add(self, reason: str, example: str = "") -> None:
        self.count[reason] += 1
        self.example.setdefault(reason, example)

    def lines(self) -> list[str]:
        return [f"{self.step}: {reason}: {n}" + (f" (first: {self.example[reason]})" if self.example[reason] else "")
                for reason, n in sorted(self.count.items())]

    def emit(self) -> None:
        for line in self.lines():
            print(line, file=sys.stderr)


def write_tsv(df: pd.DataFrame, path: Path, keys: list[str], ascending: list[bool] | bool = True) -> Path:
    """The table sorted by its key columns (a stable sort, so ties keep no accidental order: the keys must be unique),
    floats with six significant digits, missing values as empty cells."""
    if keys and df.duplicated(keys).any():
        raise ValueError(f"{path.name}: the sort keys {keys} do not identify the rows")
    out = df.sort_values(keys, ascending=ascending, kind="mergesort") if keys else df
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False, float_format="%.6g", na_rep="", lineterminator="\n")
    return path


def sign(x) -> str:
    """+, - or 0 of a number; unknown when there is none."""
    if x is None or pd.isna(x):
        return UNKNOWN
    return "+" if x > 0 else "-" if x < 0 else "0"
