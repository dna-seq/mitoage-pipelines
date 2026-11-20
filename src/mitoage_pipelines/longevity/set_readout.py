"""H27's readout: every result of the fast analyses on a restricted species set, beside the same result on all
species, with the reading fixed in H27's entry.

For a result that passed its threshold (q < 0.05) on the unrestricted run:
    reverses             the sign differs and the threshold is passed on the set
    holds                the original sign, and the threshold is passed
    holds, less precise  the original sign, the threshold not passed, the 95% interval contains the original estimate
    moves                the threshold not passed and the 95% interval excludes the original estimate
    unclear              what is left: the other sign, not passed, the interval contains the original estimate
    not assessed         the set gives no such row, or no estimate
The entry's words leave their order open where two apply (a result that passes and whose interval excludes the
original): "reverses" and "holds" are read first, and the column `interval_excludes_original` says the other.
A result that did not pass on all species and passes on a set is "new on the restricted set", never claimed.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ALPHA = 0.05
Z = 1.959964
# table (under output/longevity) -> (key columns, estimate, standard error or None when it is estimate / t)
TABLES: dict[str, tuple[list[str], str, str | None]] = {
    **{f"aa_pgls.{c}.tsv": (["class", "method", "model", "aa"], "slope", "se") for c in ("Mammalia", "Aves", "Actinopterygii", "Reptilia")},
    **{f"pgls.{c}.tsv": (["gene", "column", "residue"], "effect_log10_years", None) for c in ("Mammalia", "Aves")},
    "rates.tsv": (["class", "gene", "normalisation", "includes_saturated"], "slope", "se"),
    "motif_pgls.tsv": (["class", "subset", "feature"], "slope_per_sd", "se"),
    "generation_time.female.tsv": (["class", "maturity", "competitor", "tree", "model", "term", "feature"], "slope_sd", "se"),
    "generation_time.any.tsv": (["class", "maturity", "competitor", "tree", "model", "term", "feature"], "slope_sd", "se"),
    "thermal/composition.tsv": (["class", "model", "term", "feature"], "slope_sd", "se"),
    "thermal/pair_fits.tsv": (["class", "model", "pairing", "response", "term"], "slope", "se"),
}


def reading(est0: float, q0: float, est: float, se: float, q: float) -> str:
    """The reading of one result on a set, from its unrestricted estimate and q and the set's estimate, se and q."""
    if not np.isfinite(est) or not np.isfinite(q):
        return "not assessed"
    passed0, passed = np.isfinite(q0) and q0 < ALPHA, q < ALPHA
    if not passed0:
        return "new on the restricted set" if passed else "not passed on either"
    same = np.sign(est) == np.sign(est0)
    if passed:
        return "holds" if same else "reverses"
    if not np.isfinite(se):
        return "not assessed"
    inside = abs(est - est0) <= Z * se
    if same:
        return "holds, less precise" if inside else "moves"
    return "unclear" if inside else "moves"


def load(root: Path, table: str) -> pd.DataFrame | None:
    """One result table of a run as key columns, n, est, se, q; None when the run has no such table."""
    f = root / "output" / "longevity" / table
    if not f.is_file():
        return None
    keys, est, se = TABLES[table]
    d = pd.read_csv(f, sep="\t")
    out = d[keys].astype(str)
    out["n"] = d["n"] if "n" in d.columns else np.nan
    out["est"] = pd.to_numeric(d[est], errors="coerce")
    out["se"] = pd.to_numeric(d[se], errors="coerce") if se else (out.est / pd.to_numeric(d["t"], errors="coerce")).abs()
    out["q"] = pd.to_numeric(d["q"], errors="coerce")
    if out.duplicated(keys).any():
        raise SystemExit(f"{f}: its key columns {keys} do not name one row each")
    return out


def compare(base: Path, sets: list[str]) -> pd.DataFrame:
    """Every row of every table of `base`/all beside the same row of each set, with its reading."""
    rows = []
    for table, (keys, _, _) in TABLES.items():
        a = load(base / "all", table)
        if a is None:
            raise SystemExit(f"the unrestricted run has no {table}")
        for s in sets:
            b = load(base / s, table)
            m = a.merge(b if b is not None else a.iloc[0:0], on=keys, how="left", suffixes=("_all", ""))
            m["reading"] = [reading(*t) for t in zip(m.est_all, m.q_all, m.est, m.se, m.q)]
            m["interval_excludes_original"] = (m.est - m.est_all).abs() > Z * m.se
            m.insert(0, "result", m[keys].agg(" | ".join, axis=1))
            m.insert(0, "set", s)
            m.insert(0, "table", table)
            rows.append(m[["table", "set", "result", "n_all", "est_all", "q_all", "n", "est", "se", "q", "reading", "interval_excludes_original"]])
    return pd.concat(rows, ignore_index=True)


def main(base: Path, out: Path, sets: list[str]) -> pd.DataFrame:
    """Writes `out`/readout.tsv (every row) and `out`/summary.tsv (readings counted per table and set)."""
    out.mkdir(parents=True, exist_ok=True)
    d = compare(base, sets)
    d.to_csv(out / "readout.tsv", sep="\t", index=False, float_format="%.6g")
    s = d.groupby(["table", "set", "reading"]).size().unstack(fill_value=0).reset_index()
    s.to_csv(out / "summary.tsv", sep="\t", index=False)
    return d
