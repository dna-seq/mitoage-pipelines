"""Positive controls: amino-acid content of the 13 mt proteins against maximum lifespan.

Published effects to recover: fewer cysteines (Moosmann & Behl 2008), more threonine (Kitazoe 2008),
methionine (Aledo 2011) in long-lived animals. Relatedness is handled three ways, none a tree yet:
  raw      plain regression on log10 lifespan
  within   both sides centred on the family mean (only families with >= 3 species)
  mixed    random intercepts for order and family
each with and without log10 adult body mass.

Writes composition.tsv (one row per species: the table every later step reads) and controls.tsv.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import statsmodels.formula.api as smf
from scipy import stats

from ..genes import MT_GENES
from ..seqio import read_fasta

AAS = "ACDEFGHIKLMNPQRSTVWY"


def fit(d, aa, mass):
    """-> dict of slope per log10 lifespan (percentage points), p, n for the three treatments"""
    d = d.dropna(subset=["lm"]) if mass else d
    out = {"n": len(d)}
    if len(d) < 30: return None
    f = f"{aa} ~ ll" + (" + lm" if mass else "")
    m = smf.ols(f, d).fit(); out["raw"] = (m.params["ll"], m.pvalues["ll"])
    big = d.groupby("family").filter(lambda x: len(x) >= 3)
    if big.family.nunique() >= 5:
        c = big.copy()
        for col in [aa, "ll", "lm"] if mass else [aa, "ll"]:
            c[col] = c[col] - c.groupby("family")[col].transform("mean")
        m = smf.ols(f + " - 1", c).fit()
        # degrees of freedom lost to the family means
        df = len(c) - c.family.nunique() - (2 if mass else 1)
        t = m.params["ll"] / (m.bse["ll"] * np.sqrt(m.df_resid / df))
        out["within"] = (m.params["ll"], 2 * stats.t.sf(abs(t), df)); out["n_within"] = len(c)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = smf.mixedlm(f, d, groups="order", re_formula="1", vc_formula={"family": "0 + C(family)"}).fit(reml=True)
        out["mixed"] = (m.params["ll"], m.pvalues["ll"])
    except Exception as e:
        out["mixed_error"] = type(e).__name__
    return out


def composition(data: Path) -> pd.DataFrame:
    """Per species with all 13 proteins and a lifespan: percentage of each amino acid, log10 lifespan and mass."""
    counts, ngenes = {}, {}
    for g in MT_GENES:
        for slug, seq in read_fasta(data / "alignments" / f"{g}.aa.fasta.gz").items():
            seq = seq.replace("-", "")
            c = counts.setdefault(slug, np.zeros(len(AAS)))
            for i, a in enumerate(AAS): c[i] += seq.count(a)
            ngenes[slug] = ngenes.get(slug, 0) + 1
    sp = pl.read_parquet(data / "species.parquet").select(
        "slug", "class", "order", "family", "max_longevity_years", "adult_weight_g").to_pandas()
    sp = sp[sp.slug.map(lambda s: ngenes.get(s, 0) == 13) & sp.max_longevity_years.gt(0)].copy()
    for i, a in enumerate(AAS):
        sp["aa_" + a] = sp.slug.map(lambda s: 100 * counts[s][i] / counts[s].sum())
    sp["ll"] = np.log10(sp.max_longevity_years)
    sp["lm"] = np.log10(sp.adult_weight_g.where(sp.adult_weight_g > 0))
    return sp


def main(data: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    sp = composition(data)
    sp.to_csv(out / "composition.tsv", sep="\t", index=False)
    print(f"species with all 13 proteins and a lifespan: {len(sp)}; with body mass: {sp.lm.notna().sum()}")

    rows = []
    groups = [("all", sp)] + [(c, sp[sp["class"] == c]) for c in ["Mammalia", "Aves", "Actinopterygii", "Reptilia"]]
    for gname, d in groups:
        for aa in AAS:
            for mass in (False, True):
                r = fit(d, "aa_" + aa, mass)
                if r is None: continue
                row = {"group": gname, "aa": aa, "mass": int(mass), "n": r["n"], "n_within": r.get("n_within", "")}
                for k in ("raw", "within", "mixed"):
                    row[f"{k}_slope"], row[f"{k}_p"] = r.get(k, (np.nan, np.nan))
                row["mixed_error"] = r.get("mixed_error", "")
                rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(out / "controls.tsv", sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 250, "display.max_rows", 500):
        print(res[res.aa.isin(list("CMT"))].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
        print("\n# all 20, mammals, with mass, sorted by the within-family p")
        print(res[(res.group == "Mammalia") & (res.mass == 1)].sort_values("within_p").to_string(index=False, float_format=lambda x: f"{x:.3g}"))
