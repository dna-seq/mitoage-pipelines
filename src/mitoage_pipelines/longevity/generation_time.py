"""H11: is the mtDNA composition signal that survives the controls about lifespan, or about generation time?

The base composition of mtDNA, its GC content and the T-to-C shift at third codon positions, is the one
family of features whose correlation with lifespan survives phylogenetic contrasts and body mass in the
site's correlation table, and in mammals only once the contrast test is calibrated (the site repository's
docs/PIC_CALIBRATION.md). The A>G mutational signature on the heavy strand that would produce such a shift
has been tied to generation time across mammals (Mikhailova et al. 2022, NAR 50:10264). Lifespan and
generation time are strongly correlated, so a correlation with one can be the other's. This test puts both
in one model.

Generation time is not in the build; AnAge's age at female maturity is, and it is the proxy used here
(`--maturity any` takes the male age where the female one is missing). Per feature, standardised so that a
slope is in standard deviations of the feature, with the same species in every model:
  lifespan+mass   feature ~ log10 lifespan + log10 adult mass                     term: lifespan
  maturity+mass   feature ~ log10 female maturity + log10 adult mass              term: maturity
  both            feature ~ log10 lifespan + log10 female maturity + log10 mass   terms: lifespan, maturity
Each is a phylogenetic regression on the dated tree (pgls.py), with Pagel's lambda estimated per feature and
model ("pagel") and fixed at 1 ("brownian", the model under which the regression equals the one on
phylogenetically independent contrasts, so it can be compared with the site's contrast test). Fourfold
composition is not a covariate here, unlike H1 and H8: composition is what is being explained.

The features were fixed before any result was seen: GC content, GC skew and AT skew of the protein-coding
genes and of the whole record; A, C, G and T at fourfold-degenerate sites (`fourfold`); the share of
pyrimidine-ending codons that end in C over the protein-coding genes (`c3_share`); and the 32 codons ending
in U or C, each as a percentage of all complete codons, as the site reports them. q: Benjamini-Hochberg over
these features, within class, tree, model and term.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import linalg, optimize, stats
from statsmodels.stats.multitest import multipletests

from .aa_pgls import class_covariance, default_dataset
from .pgls import V, gls, nll

MIN_SPECIES = 25
BASES = "UCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
PYRIMIDINE_ENDING = [c for c in CODONS if c[2] in "UC"]
MODELS = {"lifespan+mass": ["ll", "lm"], "maturity+mass": ["lmat", "lm"], "both": ["ll", "lmat", "lm"]}
TERMS = {"ll": "lifespan", "lmat": "maturity"}


def traits(species: pd.DataFrame, maturity: str = "female") -> pd.DataFrame:
    """log10 lifespan, mass and age at maturity (days). `maturity` 'any': male age where female is missing."""
    log = lambda s: np.log10(s.where(s > 0))
    mat = species.female_maturity_days
    if maturity == "any":
        mat = mat.fillna(species.male_maturity_days)
    elif maturity != "female":
        raise ValueError(f"maturity must be 'female' or 'any', not {maturity!r}")
    return pd.DataFrame({"slug": species.slug, "ll": log(species.max_longevity_years),
                         "lm": log(species.adult_weight_g), "lmat": log(mat)})


def features(data: Path, fourfold: Path) -> pd.DataFrame:
    """The fixed feature set, one row per species (slug)."""
    rs = pd.read_parquet(data / "region_stats.parquet")
    rs = rs[rs.record_status.ne("retired") & rs.region.isin(["total_pc", "total_mtDNA"])]
    wide = rs.pivot_table(index="slug", columns="region", values=["gc_percent", "gc_skew", "at_skew"], aggfunc="first")
    wide.columns = [f"{v}:{r}" for v, r in wide.columns]
    cu = pd.read_parquet(data / "codon_usage.parquet")
    cu = cu[cu.region.eq("total_pc")].set_index("slug")
    total = cu[CODONS].sum(axis=1)
    cod = cu[PYRIMIDINE_ENDING].div(total, axis=0).mul(100).add_prefix("codon:")
    ending_c = cu[[c for c in CODONS if c[2] == "C"]].sum(axis=1)
    ending_u = cu[[c for c in CODONS if c[2] == "U"]].sum(axis=1)
    cod["c3_share"] = ending_c / (ending_c + ending_u)
    ff = pd.read_csv(fourfold, sep="\t").set_index("slug")[["A4", "C4", "T4"]]
    ff["G4"] = 1 - ff.A4 - ff.C4 - ff.T4
    return wide.join(cod, how="outer").join(ff, how="outer").reset_index()


def fit_terms(C: np.ndarray, y: np.ndarray, X: np.ndarray, names: list[str], lam: float | None) -> dict:
    """GLS of y on an intercept and X under Pagel's lambda (`lam` None: maximum likelihood in [0, 1]).
    -> lambda, and slope, se, t, p of every column of X, keyed by its name."""
    X1 = np.column_stack([np.ones(len(y)), X])
    if lam is None:
        lam = float(optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(X1, y, C), method="bounded").x)
    b, se, df = gls(X1, y, C, lam)
    out = {"lambda": lam, "df": df}
    for j, n in enumerate(names, start=1):
        t = b[j] / se[j]
        out[n] = (b[j], se[j], t, 2 * stats.t.sf(abs(t), df))
    return out


def whitened_vif(C: np.ndarray, X: np.ndarray, lam: float = 1.0) -> np.ndarray:
    """Variance inflation of each column of X in the GLS design (intercept included, then dropped)."""
    L = linalg.cholesky(V(lam, C), lower=True)
    Xw = linalg.solve_triangular(L, np.column_stack([np.ones(len(X)), X]), lower=True)
    inv = np.linalg.inv(Xw.T @ Xw)
    return np.array([inv[j, j] * (Xw[:, j] @ Xw[:, j] - (Xw[:, j] @ Xw[:, 0]) ** 2 / (Xw[:, 0] @ Xw[:, 0]))
                     for j in range(1, Xw.shape[1])])


def select(data: Path, fourfold: Path, tree_file: Path, cls: str, maturity: str) -> tuple[pd.DataFrame, list[tuple[str, int]]]:
    sp = pd.read_parquet(data / "species.parquet")
    sp = sp[sp["class"] == cls]
    steps = [(f"{cls} in the build", len(sp))]
    sp = sp[default_dataset(sp)]
    steps.append(("default dataset (live; lifespan quality high or acceptable; not estimated)", len(sp)))
    tr = traits(sp, maturity).merge(features(data, fourfold), on="slug", how="left")
    tr = tr.dropna(subset=["ll", "lm"])
    steps.append(("and body mass", len(tr)))
    tr = tr.dropna(subset=["lmat"])
    steps.append((f"and age at maturity ({maturity})", len(tr)))
    present = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    tr = tr[tr.slug.isin(present)].reset_index(drop=True)
    steps.append(("and a tip in the dated tree", len(tr)))
    return tr, steps


def run_class(data: Path, fourfold: Path, tree_file: Path, cls: str, maturity: str) -> tuple[list[dict], dict]:
    sp, steps = select(data, fourfold, tree_file, cls, maturity)
    for label, n in steps:
        print(f"{cls}: {label}: {n}")
    if len(sp) < MIN_SPECIES:
        print(f"{cls}: fewer than {MIN_SPECIES} species, not tested")
        return [], {}
    C, order, extra = class_covariance(tree_file, list(sp.slug))
    sp = sp.set_index("slug").loc[order].reset_index()
    vif = whitened_vif(C, sp[["ll", "lmat", "lm"]].to_numpy())
    info = {"class": cls, "n": len(sp), "r_lifespan_maturity": float(np.corrcoef(sp.ll, sp.lmat)[0, 1]),
            "vif_lifespan": vif[0], "vif_maturity": vif[1], "vif_mass": vif[2]}
    print(f"{cls}: {len(sp)} species ({extra} other tips under their ancestor ignored); r(lifespan, maturity) across "
          f"species {info['r_lifespan_maturity']:.2f}; VIF on the tree (lambda 1): lifespan {vif[0]:.1f}, "
          f"maturity {vif[1]:.1f}, mass {vif[2]:.1f}")
    feats = [c for c in sp.columns if c not in ("slug", "ll", "lm", "lmat")]
    rows = []
    for f in feats:
        x = sp[f].to_numpy(dtype=float)
        ok = np.isfinite(x)
        if ok.sum() < MIN_SPECIES or np.std(x[ok]) == 0:
            continue
        y = (x[ok] - x[ok].mean()) / x[ok].std()
        Cs = C[np.ix_(ok, ok)]
        for tree, lam in (("pagel", None), ("brownian", 1.0)):
            for model, cols in MODELS.items():
                r = fit_terms(Cs, y, sp.loc[ok, cols].to_numpy(), cols, lam)
                for col in cols:
                    if col not in TERMS:
                        continue
                    b, se, t, p = r[col]
                    rows.append({"class": cls, "maturity": maturity, "tree": tree, "model": model, "term": TERMS[col],
                                 "feature": f, "n": int(ok.sum()), "lambda": r["lambda"], "slope_sd": b, "se": se,
                                 "t": t, "df": r["df"], "p": p})
    return rows, info


def main(data: Path, fourfold: Path, tree_file: Path, out: Path, classes: tuple[str, ...] = ("Mammalia", "Aves"),
         maturity: str = "female") -> pd.DataFrame:
    """Writes generation_time.<maturity>.tsv (every fit) and generation_time.<maturity>.classes.tsv
    (per class: species, the correlation of lifespan and maturity, and their variance inflation)."""
    out.mkdir(parents=True, exist_ok=True)
    rows, infos = [], []
    for cls in classes:
        r, info = run_class(data, fourfold, tree_file, cls, maturity)
        rows += r
        if info:
            infos.append(info)
    if not rows:
        raise SystemExit("no class had enough species with lifespan, body mass, maturity and a tip in the tree")
    res = pd.DataFrame(rows)
    res["q"] = res.groupby(["class", "tree", "model", "term"]).p.transform(
        lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index))
    res = res.sort_values(["class", "tree", "model", "term", "p"]).reset_index(drop=True)
    res.to_csv(out / f"generation_time.{maturity}.tsv", sep="\t", index=False, float_format="%.5g")
    pd.DataFrame(infos).to_csv(out / f"generation_time.{maturity}.classes.tsv", sep="\t", index=False, float_format="%.4g")
    summary = (res[res.tree.eq("pagel")].assign(sig=lambda d: d.q < 0.05)
               .groupby(["class", "model", "term"]).agg(tests=("p", "size"), q_below_005=("sig", "sum"),
                                                        smallest_p=("p", "min")).reset_index())
    with pd.option_context("display.width", 200):
        print(summary.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res
