"""H31: the table of a masked DNA model's likelihoods along one circular genome (numpy only; the model is a callable).

Every place is hidden in turn, in groups of places `stride` apart, and the model's log-probabilities of A, C, G and
T at the hidden places are kept. The genome is read twice, in the record's order and from its middle, and each place
takes its value from the reading in which it lies farther from the ends of the sequence.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

BASES = "ACGT"
N, MASK = 4, 5
STRIDE = 128
Predict = Callable[[np.ndarray], np.ndarray]     # codes [B, L] with MASK at hidden places -> log-probabilities [B, L, 4]


def encode(seq: str) -> np.ndarray:
    """A 0, C 1, G 2, T 3, anything else 4 (unknown)."""
    lut = np.full(256, N, dtype=np.int64)
    for i, b in enumerate(BASES):
        lut[ord(b)] = i; lut[ord(b.lower())] = i
    return lut[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def padded(codes: np.ndarray, multiple: int) -> np.ndarray:
    """The codes followed by unknowns up to the next multiple."""
    out = np.full(-(-len(codes) // multiple) * multiple, N, dtype=np.int64)
    out[:len(codes)] = codes
    return out


def groups(n: int, stride: int = STRIDE) -> list[np.ndarray]:
    """The places 0..n-1 in sets whose members are `stride` apart; every place is in exactly one."""
    return [np.arange(k, n, stride) for k in range(min(stride, n))]


def one_reading(codes: np.ndarray, predict: Predict, stride: int, multiple: int, batch: int) -> np.ndarray:
    """[n, 4] log-probabilities, each place hidden once; a place that holds an unknown is not scored (NaN)."""
    n = len(codes)
    out = np.full((n, 4), np.nan)
    sets = groups(n, stride)
    for a in range(0, len(sets), batch):
        part = sets[a:a + batch]
        rows = np.stack([padded(codes, multiple) for _ in part])
        for r, places in zip(rows, part):
            r[places] = MASK
        got = predict(rows)
        if got.shape != (len(part), rows.shape[1], 4):
            raise ValueError(f"the model returned {got.shape} for {rows.shape} codes")
        for g, places in zip(got, part):
            out[places] = g[places]
    out[codes == N] = np.nan
    return out


def table(codes: np.ndarray, predict: Predict, stride: int = STRIDE, multiple: int = 128, batch: int = 8) -> np.ndarray:
    """[n, 4] log-probabilities of A, C, G, T at each place of a circular genome, from two readings (the record's
    order, and from the middle), each place from the reading where it lies farther from the nearer end."""
    n = len(codes)
    half = n // 2
    first = one_reading(codes, predict, stride, multiple, batch)
    turned = one_reading(np.roll(codes, -half), predict, stride, multiple, batch)
    second = np.roll(turned, half, axis=0)                    # back to the record's places
    i = np.arange(n)
    edge_first = np.minimum(i, n - 1 - i)
    j = (i - half) % n                                        # the place's index in the turned reading
    edge_second = np.minimum(j, n - 1 - j)
    return np.where((edge_first >= edge_second)[:, None], first, second)


def write(path: Path, seq: str, tab: np.ndarray) -> None:
    """pos (1-based), ref, lp_A, lp_C, lp_G, lp_T; an unscored place has empty values."""
    if len(seq) != len(tab):
        raise ValueError("the table and the sequence differ in length")
    with open(path, "w") as f:
        f.write("pos\tref\t" + "\t".join(f"lp_{b}" for b in BASES) + "\n")
        for k, (b, row) in enumerate(zip(seq.upper(), tab), 1):
            f.write(f"{k}\t{b}\t" + "\t".join("" if np.isnan(v) else f"{v:.5f}" for v in row) + "\n")
