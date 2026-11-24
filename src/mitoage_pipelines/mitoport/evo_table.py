"""H35, part 1: the table of a left-to-right DNA model's likelihoods along one circular genome (numpy only; the
model is a callable).

A left-to-right model gives each base from the bases before it. Each strand is read twice, from the record's start
and from its middle, and a place takes its value from the reading in which more bases precede it (at least half the
genome). The minus strand is read as the reverse complement and turned back: the plus strand's base b at place p is
the complement of b at place n - 1 - p there. The two strands' log-probabilities are averaged and brought back to
sum to one over A, C, G and T; a place that holds another letter is not scored.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

try:
    from .kron_table import BASES, N, encode
except ImportError:          # read as one file by slurm/h35-evo2-table.py, outside the package
    from kron_table import BASES, N, encode

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
Predict = Callable[[str], np.ndarray]     # a sequence -> [n, 4] log-probabilities of A, C, G, T given the bases before; row 0 is not used


def one_strand(seq: str, predict: Predict) -> np.ndarray:
    """[n, 4]: each place read with the bases before it on this strand, from the start or from the middle of the
    circle, whichever puts more bases before it."""
    n, half = len(seq), len(seq) // 2
    first = np.asarray(predict(seq), dtype=float)
    turned = np.asarray(predict(seq[half:] + seq[:half]), dtype=float)
    for got in (first, turned):
        if got.shape != (n, 4):
            raise ValueError(f"the model returned {got.shape} for {n} bases")
    second = np.roll(turned, half, axis=0)                 # row p of `second` is place p, read with p - half (mod n) bases before
    before_first = np.arange(n)
    before_second = (np.arange(n) - half) % n
    return np.where((before_first >= before_second)[:, None], first, second)


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def turned_back(minus: np.ndarray) -> np.ndarray:
    """The table of the reverse complement, as the plus strand's places and bases: rows reversed, A with T and C with
    G exchanged."""
    return minus[::-1, ::-1]


def normalised(lp: np.ndarray) -> np.ndarray:
    m = lp.max(axis=1, keepdims=True)
    return lp - (m + np.log(np.exp(lp - m).sum(axis=1, keepdims=True)))


def tables(seq: str, predict: Predict) -> dict[str, np.ndarray]:
    """{"plus": the plus strand's reading, "minus": the minus strand's turned back, "both": their mean, normalised};
    each [n, 4], NaN where the genome holds no A, C, G or T."""
    seq = seq.upper()
    plus = normalised(one_strand(seq, predict))
    minus = normalised(turned_back(one_strand(reverse_complement(seq), predict)))
    out = {"plus": plus, "minus": minus, "both": normalised((plus + minus) / 2)}
    unknown = encode(seq) == N
    for t in out.values():
        t[unknown] = np.nan
    return out


__all__ = ["BASES", "tables", "one_strand", "reverse_complement", "turned_back", "normalised"]
