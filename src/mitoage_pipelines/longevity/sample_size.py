"""H27, sample size as a covariate: the lifespan slopes of H1 and H11 with and without AnAge's sample size in the
model, on the same species.

A maximum lifespan rises with the number of animals it was read from, so a slope on lifespan could carry how
well a species is sampled. AnAge grades its sample in five ordered classes; here it enters each model as one
ordered term (tiny 0, small 1, medium 2, large 3, huge 4), on the species of the set `not_unreliable` whose
lifespan is AnAge's and has a sample size. Each model is fitted twice on those species, without the term and
with it, so the two slopes differ by the term alone.

  H1    aa ~ log10 lifespan + log10 mass [+ A4 + C4 + T4], Pagel's lambda, and the mixed model (aa_pgls.py)
  H11   feature (standardised) ~ the entry's models that hold a lifespan term (generation_time.py); the two
        slow-fast axes are not refitted: lifespan is inside the axis, not a term beside it

q: Benjamini-Hochberg within the groups of the entries (H1: class, method, model; H11: class, competitor, tree,
model), separately without and with the term. The reading of a result that passes without the term is
`set_readout.reading`'s, with the fit without the term as the original.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from . import aa_pgls, generation_time as gt, species_sets
from .set_readout import ALPHA, reading

SET = "not_unreliable"
ORDER = {"tiny": 0, "small": 1, "medium": 2, "large": 3, "huge": 4}
KEYS = ["analysis", "class", "method", "competitor", "model", "feature"]


def sample_size(species: pd.DataFrame) -> pd.DataFrame:
    """slug and ss (0 to 4) of the species whose lifespan is AnAge's with a sample size. A sample size that is
    none of the five classes is refused: it would otherwise drop out as if the species had none."""
    for c in ("slug", "lifespan_source", "lifespan_sample_size"):
        if c not in species.columns:
            raise SystemExit(f"the species table lacks {c}")
    own = species[species.lifespan_source.eq("AnAge") & species.lifespan_sample_size.notna()]
    unknown = sorted(set(own.lifespan_sample_size) - set(ORDER))
    if unknown:
        raise SystemExit(f"lifespan_sample_size holds values that are no sample size: {', '.join(map(str, unknown))}")
    return pd.DataFrame({"slug": own.slug, "ss": own.lifespan_sample_size.map(ORDER).astype(float)}).reset_index(drop=True)


def both(C: np.ndarray, y: np.ndarray, X: np.ndarray, ss: np.ndarray, lam: float | None = None) -> dict:
    """The slope of the first column of X without and with `ss` as one more column, on the same rows."""
    names = [f"x{j}" for j in range(X.shape[1])]
    a = gt.fit_terms(C, y, X, names, lam)
    b = gt.fit_terms(C, y, np.column_stack([X, ss]), [*names, "ss"], lam)
    out = {}
    for tag, r in (("without", a), ("with", b)):
        slope, se, _, p = r["x0"]
        out |= {f"slope_{tag}": slope, f"se_{tag}": se, f"p_{tag}": p, f"lambda_{tag}": r["lambda"]}
    out["ss_slope"], out["ss_se"], _, out["ss_p"] = b["ss"]
    return out


def h1(data: Path, composition: Path, fourfold: Path, tree_file: Path, ss: pd.DataFrame, cls: str) -> list[dict]:
    sp, _ = aa_pgls.select(data, composition, fourfold, tree_file, cls)
    n0 = len(sp)
    sp = sp.merge(ss, on="slug")
    C, order, _ = aa_pgls.class_covariance(tree_file, list(sp.slug))
    sp = sp.set_index("slug").loc[order].reset_index()
    print(f"H1, {cls}: {n0} species in the set's test, {len(sp)} of them with an AnAge sample size "
          + ", ".join(f"{k} {int((sp.ss == v).sum())}" for k, v in ORDER.items()))
    rows = []
    for model, covs in aa_pgls.MODELS.items():
        X = np.column_stack([sp.ll.to_numpy(), sp[covs].to_numpy()])
        for aa in aa_pgls.AAS:
            base = {"analysis": "H1", "class": cls, "competitor": "", "model": model, "feature": aa}
            rows.append({**base, "method": "pgls", "n": len(sp), **both(C, sp["aa_" + aa].to_numpy(), X, sp.ss.to_numpy())})
            a, b = aa_pgls.fit_mixed(sp, aa, covs), aa_pgls.fit_mixed(sp, aa, [*covs, "ss"])
            rows.append({**base, "method": "mixed", "n": a["n"], "slope_without": a["slope"], "se_without": a["se"], "p_without": a["p"],
                         "slope_with": b["slope"], "se_with": b["se"], "p_with": b["p"]})
    return rows


def h11_models() -> dict[str, list[tuple[str, list[str]]]]:
    """competitor -> the entry's models that hold a lifespan term."""
    out: dict[str, list[tuple[str, list[str]]]] = {}
    for comp, model, cols in gt.models():
        if "ll" in cols:
            out.setdefault(comp, []).append((model, cols))
    return out


def h11(data: Path, fourfold: Path, tree_file: Path, ss: pd.DataFrame, cls: str, maturity: str) -> list[dict]:
    sp, _ = gt.select(data, fourfold, tree_file, cls, maturity)
    sp = sp.merge(ss, on="slug")
    rows = []
    for comp, specs in h11_models().items():
        cols = sorted({c for _, cs in specs for c in cs}, key=gt.TRAIT_COLUMNS.index)
        sub = sp.dropna(subset=cols).reset_index(drop=True)
        if len(sub) < gt.MIN_SPECIES:
            print(f"H11, {cls}, {comp}: {len(sub)} species with an AnAge sample size, fewer than {gt.MIN_SPECIES}: not tested")
            continue
        C, order, _ = aa_pgls.class_covariance(tree_file, list(sub.slug))
        sub = sub.set_index("slug").loc[order].reset_index()
        print(f"H11, {cls}, {comp}: {len(sub)} species with an AnAge sample size")
        feats = [c for c in sub.columns if c not in ["slug", "ss", *gt.TRAIT_COLUMNS]]
        for f in feats:
            x = sub[f].to_numpy(dtype=float)
            ok = np.isfinite(x)
            if ok.sum() < gt.MIN_SPECIES or np.std(x[ok]) == 0:
                continue
            y = (x[ok] - x[ok].mean()) / x[ok].std()
            Cs = C[np.ix_(ok, ok)]
            for tree, lam in (("pagel", None), ("brownian", 1.0)):
                for model, mcols in specs:
                    rows.append({"analysis": "H11", "class": cls, "method": tree, "competitor": comp, "model": model, "feature": f,
                                 "n": int(ok.sum()), **both(Cs, y, sub.loc[ok, mcols].to_numpy(), sub.ss.to_numpy()[ok], lam)})
    return rows


def read(res: pd.DataFrame) -> pd.DataFrame:
    """q without and with the term within the entries' groups, and the reading of each row."""
    res = res.copy()
    groups = ["analysis", "class", "method", "competitor", "model"]
    for tag in ("without", "with"):
        res[f"q_{tag}"] = res.groupby(groups)[f"p_{tag}"].transform(
            lambda p: pd.Series(multipletests(p.dropna(), method="fdr_bh")[1], index=p.dropna().index).reindex(p.index))
    res["reading"] = [reading(r.slope_without, r.q_without, r.slope_with, r.se_with, r.q_with) for r in res.itertuples()]
    return res


def summary(res: pd.DataFrame) -> pd.DataFrame:
    """Per analysis, class, method and model: tests, how many pass without the term, their readings with it, the
    median of slope with / slope without over those, and how often the term itself has p < 0.05."""
    rows = []
    for key, d in res.groupby(["analysis", "class", "method", "competitor", "model"], sort=False):
        passed = d[d.q_without < ALPHA]
        ratio = (passed.slope_with / passed.slope_without).median() if len(passed) else np.nan
        counts = passed.reading.value_counts()
        rows.append({**dict(zip(["analysis", "class", "method", "competitor", "model"], key)), "n": int(d.n.max()), "tests": len(d),
                     "pass_without": len(passed), "holds": int(counts.get("holds", 0)),
                     "holds_less_precise": int(counts.get("holds, less precise", 0)), "moves": int(counts.get("moves", 0)),
                     "reverses": int(counts.get("reverses", 0)), "other": int(len(passed) - sum(counts.get(k, 0) for k in ("holds", "holds, less precise", "moves", "reverses"))),
                     "new_with": int((d.reading == "new on the restricted set").sum()),
                     "median_slope_ratio": ratio, "term_p_below_005": int((d.ss_p < ALPHA).sum()) if "ss_p" in d else 0})
    return pd.DataFrame(rows)


def main(data: Path, composition: Path, fourfold: Path, tree_file: Path, out: Path,
         classes: tuple[str, ...] = ("Mammalia", "Aves"), maturity: str = "female") -> pd.DataFrame:
    """Writes sample_size_covariate.tsv (every fit) and sample_size_covariate.summary.tsv under `out`."""
    if species_sets.active() != SET:
        raise SystemExit(f"this model is defined on the species set {SET}: set {species_sets.ENV}={SET}")
    ss = sample_size(pd.read_parquet(data / "species.parquet"))
    rows = []
    for cls in classes:
        rows += h1(data, composition, fourfold, tree_file, ss, cls)
        rows += h11(data, fourfold, tree_file, ss, cls, maturity)
    res = read(pd.DataFrame(rows))
    res = res.sort_values(KEYS).reset_index(drop=True)
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "sample_size_covariate.tsv", sep="\t", index=False, float_format="%.5g")
    s = summary(res)
    s.to_csv(out / "sample_size_covariate.summary.tsv", sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(s.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res
