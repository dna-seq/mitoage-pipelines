"""H2's power: the smallest effect of a residue on log10 lifespan that the column test detects with 80 %
power at the significance a lone discovery needs to pass 5 % false discovery (Benjamini-Hochberg: p <= 0.05 k / m
for k discoveries among m tests).

The standard error of an effect depends on how the carriers sit on the tree; it is read from the tests themselves
(effect / t in pgls.<Class>.tsv), as the median within bins of carrier count. Smallest detectable effect:
se (t[1 - alpha/2] + t[power]), on the residual degrees of freedom. A normal sampling distribution of the effect is
assumed, and lambda as estimated.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

BINS = [(15, 20), (20, 30), (30, 50), (50, 100), (100, 200), (200, 400)]


def detectable(se: float, alpha: float, df: int, power: float = 0.8) -> float:
    """The effect at which a two-sided t test at `alpha` rejects with probability `power`."""
    return se * (stats.t.ppf(1 - alpha / 2, df) + stats.t.ppf(power, df))


def main(table: Path, out: Path, discoveries: tuple[int, ...] = (1, 2), power: float = 0.8, fdr: float = 0.05) -> pd.DataFrame:
    d = pd.read_csv(table, sep="\t")
    d["se"] = d.effect_log10_years / d.t
    m = len(d); n = int(np.median(d.carriers + d.others)); df = n - 6     # intercept, mass, A4, C4, T4, residue
    rows = []
    for k in discoveries:
        alpha = fdr * k / m
        for lo, hi in BINS:
            g = d[(d.carriers >= lo) & (d.carriers < hi)]
            if g.empty: continue
            e = detectable(float(g.se.median()), alpha, df, power)
            rows.append({"discoveries": k, "alpha": alpha, "carriers": f"{lo}-{hi - 1}", "tests": len(g),
                         "median_se": g.se.median(), "detectable_log10_years": e, "detectable_fold": 10 ** e})
    res = pd.DataFrame(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False, float_format="%.4g")
    print(f"{m} tests, median {n} species per test; carriers quartiles {d.carriers.quantile([.25, .5, .75]).tolist()}")
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res
