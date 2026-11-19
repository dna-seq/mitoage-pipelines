"""H1 on the dated species tree: amino-acid content of the 13 mt proteins against maximum lifespan, as a
phylogenetic regression with Pagel's lambda.

Per amino acid, its percentage over the 13 proteins (composition.tsv, from `controls`) is the response:
  mass               aa ~ log10 lifespan + log10 adult mass
  mass+composition   aa ~ log10 lifespan + log10 adult mass + A4 + C4 + T4   (fourfold sites, from `fourfold`)
Lambda is estimated by maximum likelihood for each amino acid and model (the phylogenetic signal of the
residual), on the covariance of the dated species tree (pgls.py). The slope of log10 lifespan is in percentage
points per tenfold lifespan. q: Benjamini-Hochberg over the 20 amino acids, within class, method and model.
The earlier test of H1 (`controls`, `composition`: random intercepts of order and family) is fitted again on
the same species (method `mixed`), so that what the tree changes is told apart from what the species set does.

Species: the manuscript's default dataset (not retired; lifespan record of high or acceptable quality and not
estimated), with body mass, all 13 proteins, fourfold counts, and a tip in the dated tree.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import optimize, stats
from statsmodels.stats.multitest import multipletests

from .pgls import gls, nll, shared_path_covariance
from . import species_sets
from .. import paths

AAS = "ACDEFGHIKLMNPQRSTVWY"
QUALITY_OK = ("high", "acceptable")
MODELS = {"mass": ["lm"], "mass+composition": ["lm", "A4", "C4", "T4"]}


def default_dataset(species: pd.DataFrame) -> pd.Series:
    """True for the species of the manuscript's default dataset: not retired, lifespan quality high or
    acceptable, lifespan record not estimated. When MITOAGE_SPECIES_SET names a species set (H27), only the
    species of that set, by their slug, which `species` must then carry."""
    keep = (species.record_status.ne("retired") & species.lifespan_quality.isin(QUALITY_OK)
            & species.lifespan_record_type.ne("estimated"))
    name = species_sets.active()
    if name is None:
        return keep
    if "slug" not in species.columns:
        raise SystemExit(f"{species_sets.ENV}={name}: this analysis selects its species without their slug, "
                         "so the set cannot be applied to it")
    return keep & species.slug.isin(species_sets.slugs(paths.BUILD, name))


def fit(C: np.ndarray, y: np.ndarray, X: np.ndarray, lam: float | None = None) -> dict:
    """GLS of y on X with covariance C scaled by Pagel's lambda; an intercept is added and the first column
    of X is the one tested. `lam` None: estimated by maximum likelihood in [0, 1]."""
    X = np.column_stack([np.ones(len(y)), X])
    if lam is None:
        lam = float(optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(X, y, C), method="bounded").x)
    b, se, df = gls(X, y, C, lam)
    t = b[1] / se[1]
    return {"lambda": lam, "slope": b[1], "se": se[1], "t": t, "df": df, "p": 2 * stats.t.sf(abs(t), df)}


def fit_mixed(sp: pd.DataFrame, aa: str, covs: list[str]) -> dict:
    """The earlier H1 model: random intercepts of order and of family within it (as `composition`), on the
    species with both an order and a family."""
    import warnings
    import statsmodels.formula.api as smf
    sp = sp.dropna(subset=["order", "family"])
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = smf.mixedlm(f"aa_{aa} ~ ll + " + " + ".join(covs), sp, groups="order",
                            vc_formula={"family": "0 + C(family)"}).fit()
        return {"n": len(sp), "lambda": np.nan, "slope": m.params["ll"], "se": m.bse["ll"], "t": m.tvalues["ll"], "p": m.pvalues["ll"]}
    except Exception as e:                                      # reported, not hidden
        print(f"mixed model failed for {aa} ({' + '.join(covs)}): {type(e).__name__}")
        return {"n": len(sp), "lambda": np.nan, "slope": np.nan, "se": np.nan, "t": np.nan, "p": np.nan}


def class_covariance(tree_file: Path, slugs: list[str]) -> tuple[np.ndarray, list[str], int]:
    """Covariance of the species of `slugs` that are tips of the tree, taken under their common ancestor (as
    `pgls` does). -> (C, those slugs in the order of C, tips under the ancestor that are not in `slugs`)."""
    tree = Phylo.read(tree_file, "newick")
    wanted = set(slugs)
    anc = tree.common_ancestor([t for t in tree.get_terminals() if t.name in wanted])
    tips = {t.name: i for i, t in enumerate(anc.get_terminals())}
    keep = [s for s in slugs if s in tips]
    idx = np.array([tips[s] for s in keep])
    C = shared_path_covariance(anc, tips)[np.ix_(idx, idx)]
    return C, keep, len(tips) - len(keep)


def select(data: Path, composition: Path, fourfold: Path, tree_file: Path, cls: str) -> tuple[pd.DataFrame, list[tuple[str, int]]]:
    """The species of one class the test uses, and how many remain after each filter."""
    comp = pd.read_csv(composition, sep="\t")
    sp = comp[comp["class"] == cls].merge(
        pd.read_parquet(data / "species.parquet", columns=["slug", "record_status", "lifespan_quality", "lifespan_record_type"]),
        on="slug", how="left")
    steps = [("all 13 proteins and a lifespan", len(sp))]
    sp = sp[default_dataset(sp)]
    steps.append(("default dataset (live; lifespan quality high or acceptable; not estimated)", len(sp)))
    sp = sp.dropna(subset=["lm"])
    steps.append(("and body mass", len(sp)))
    sp = sp.merge(pd.read_csv(fourfold, sep="\t"), on="slug")
    steps.append(("and fourfold counts", len(sp)))
    present = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    sp = sp[sp.slug.isin(present)].reset_index(drop=True)
    steps.append(("and a tip in the dated tree", len(sp)))
    return sp, steps


def main(data: Path, composition: Path, fourfold: Path, tree_file: Path, out: Path, cls: str = "Mammalia") -> pd.DataFrame:
    """Writes `out`: class, method (pgls, mixed), model, aa, n, lambda, slope, se, t, p, q."""
    out.parent.mkdir(parents=True, exist_ok=True)
    sp, steps = select(data, composition, fourfold, tree_file, cls)
    for label, n in steps:
        print(f"{cls}: {label}: {n}")
    C, order, extra = class_covariance(tree_file, list(sp.slug))
    sp = sp.set_index("slug").loc[order].reset_index()
    print(f"{cls}: {len(sp)} species in the test; {extra} other tips under their common ancestor (ignored)")
    rows = []
    for model, covs in MODELS.items():
        X = np.column_stack([sp.ll.to_numpy(), sp[covs].to_numpy()])
        for aa in AAS:
            r = fit(C, sp["aa_" + aa].to_numpy(), X); r.pop("df")
            rows.append({"class": cls, "method": "pgls", "model": model, "aa": aa, "n": len(sp), **r})
            rows.append({"class": cls, "method": "mixed", "model": model, "aa": aa, **fit_mixed(sp, aa, covs)})
    res = pd.DataFrame(rows)
    res["q"] = res.groupby(["method", "model"]).p.transform(
        lambda p: pd.Series(multipletests(p.dropna(), method="fdr_bh")[1], index=p.dropna().index).reindex(p.index))
    res = res.sort_values(["method", "model", "aa"], ascending=[False, True, True])
    res.to_csv(out, sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 200):
        print(res[res.method == "pgls"].sort_values(["model", "p"]).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res
