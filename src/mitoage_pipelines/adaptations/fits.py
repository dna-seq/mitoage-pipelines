"""Phylogenetic regressions of H21-H25: composition against a trait (R3), lifespan against hibernation (H24's P,
H25's B1), lifespan's composition signal with and without hibernation (H24's L2), and the sister-pair rate against
lifespan inside bats (H25's B2).

Every fit is GLS on the build's dated tree with Pagel's lambda by maximum likelihood (`thermal.fit`), on the species
that have every term. Composition: H11's 43 features (`generation_time.features`), each standardised; q is
Benjamini-Hochberg over the features within one model and term.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from statsmodels.stats.multitest import multipletests

from ..longevity.aa_pgls import class_covariance
from ..longevity.thermal import disjoint_pairs, fit, pair_table

MIN_SPECIES = 20
MIN_BAT_PAIRS = 12
L2_LARGELY, L2_NOT = 0.5, 0.8          # H24's reading of the median ratio of lifespan slopes
L2_MIN_FEATURES = 5


def feature_fits(d: pd.DataFrame, feats: pd.DataFrame, terms: list[str], tree_file: Path, label: str) -> pd.DataFrame:
    """Each feature on `terms` over the species of `d` (slug and the terms) that have every term and the feature."""
    names = [c for c in feats.columns if c != "slug"]
    d = d.dropna(subset=terms).merge(feats, on="slug")
    rows = []
    if len(d) < MIN_SPECIES:
        return pd.DataFrame([{"model": label, "n": len(d), "status": f"not assessed: fewer than {MIN_SPECIES} species"}])
    C, order, _ = class_covariance(tree_file, list(d.slug))
    d = d.set_index("slug").loc[order].reset_index()
    X = d[terms].to_numpy(dtype=float)
    for f in names:
        x = d[f].to_numpy(dtype=float); ok = np.isfinite(x)
        if ok.sum() < MIN_SPECIES or np.std(x[ok]) == 0:
            rows.append({"model": label, "feature": f, "n": int(ok.sum()), "status": "not assessed: too few values"})
            continue
        y = (x[ok] - x[ok].mean()) / x[ok].std()
        r = fit(C[np.ix_(ok, ok)], y, X[ok], terms)
        for t in terms:
            b, se, tt, p = r[t]
            rows.append({"model": label, "term": t, "feature": f, "n": int(ok.sum()), "status": "assessed",
                         "lambda": r["lambda"], "slope_sd": b, "se": se, "t": tt, "p": p})
    res = pd.DataFrame(rows)
    ok = res.status.eq("assessed")
    res.loc[ok, "q"] = res[ok].groupby("term").p.transform(lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index))
    return res


def one_fit(d: pd.DataFrame, y: str, terms: list[str], tree_file: Path, label: str) -> list[dict]:
    d = d.dropna(subset=[y, *terms])
    if len(d) < MIN_SPECIES:
        return [{"model": label, "n": len(d), "status": f"not assessed: fewer than {MIN_SPECIES} species"}]
    C, order, _ = class_covariance(tree_file, list(d.slug))
    d = d.set_index("slug").loc[order].reset_index()
    r = fit(C, d[y].to_numpy(float), d[terms].to_numpy(float), terms)
    return [{"model": label, "response": y, "term": t, "n": len(d), "status": "assessed", "lambda": r["lambda"],
             "slope": r[t][0], "se": r[t][1], "t": r[t][2], "p": r[t][3]} for t in terms]


def l2(without: pd.DataFrame, with_hib: pd.DataFrame) -> dict:
    """H24's L2: lifespan's composition signal with hibernation in the model against without, on the same species."""
    a = without[without.status.eq("assessed") & without.term.eq("ll")].set_index("feature")
    b = with_hib[with_hib.status.eq("assessed") & with_hib.term.eq("ll")].set_index("feature")
    if a.empty or b.empty:
        return {"readout": "L2", "status": "not assessed: a model was not fitted"}
    if set(a.n) != set(b.n) or not a.index.equals(b.index):
        raise ValueError("L2 compares two models on different species or features")
    passing = a.index[a.q < 0.05]
    out = {"readout": "L2", "n_species": int(a.n.iloc[0]), "passing_without": int(len(passing)),
           "passing_with": int((b.q < 0.05).sum())}
    if len(passing) < L2_MIN_FEATURES:
        return {**out, "status": f"not assessed: {len(passing)} features pass without hibernation, fewer than {L2_MIN_FEATURES}"}
    ratio = float(np.median(np.abs(b.loc[passing, "slope_sd"]) / np.abs(a.loc[passing, "slope_sd"])))
    halved = out["passing_with"] <= out["passing_without"] / 2
    reading = ("largely hibernation's" if ratio < L2_LARGELY and halved else "not hibernation's" if ratio > L2_NOT else "partly")
    return {**out, "status": "assessed", "median_ratio": ratio, "count_halved": halved, "reading": reading}


def bat_pairs(sp: pd.DataFrame, dist, tree_file: Path) -> tuple[pd.DataFrame, list[dict]]:
    """H25's B2 by H20's design: the largest set of bat pairs with no shared branch; per pair, log10 dS per million
    years between the sisters on the pair's mean log10 lifespan and log10 mass (GLS, Pagel's lambda, each pair at its
    first species' tip). dN per Myr and dN/dS beside it."""
    tree = Phylo.read(tree_file, "newick")
    members = set(sp.dropna(subset=["ll", "lm"]).slug) & {t.name for t in tree.get_terminals()} & dist.aligned
    p = pair_table(disjoint_pairs(tree, members), dist, sp, ["ll", "lm"])
    d = p[p.usable].reset_index(drop=True)
    base = {"model": "B2", "n_species": len(members), "n_pairs_all": len(p), "n_saturated": int(p.saturated.sum()),
            "n_pairs": len(d)}
    if len(d) < MIN_BAT_PAIRS:
        return p, [{**base, "status": f"not assessed: {len(d)} pairs, fewer than {MIN_BAT_PAIRS}"}]
    C, order, _ = class_covariance(tree_file, list(d.sp_a))
    d = d.set_index("sp_a").loc[order].reset_index()
    rows = []
    for resp, col in (("log10_dS_per_myr", "dS_per_myr"), ("log10_dN_per_myr", "dN_per_myr"), ("log10_dNdS", "dNdS")):
        r = fit(C, np.log10(d[col].to_numpy(float)), d[["ll", "lm"]].to_numpy(float), ["ll", "lm"])
        for t in ("ll", "lm"):
            rows.append({**base, "response": resp, "term": t, "status": "assessed", "lambda": r["lambda"],
                         "slope": r[t][0], "se": r[t][1], "t": r[t][2], "df": r["df"], "p": r[t][3]})
    return p, rows
