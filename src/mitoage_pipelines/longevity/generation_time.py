"""H11: how the mtDNA composition signal that survives the controls is shared across life-history traits.

The base composition of mtDNA, its GC content and the T-to-C shift at third codon positions, is the one
family of features whose correlation with lifespan survives phylogenetic contrasts and body mass in the
site's correlation table, and in mammals only once the contrast test is calibrated (the site repository's
docs/PIC_CALIBRATION.md). Two traits that go with lifespan have a claim on it. Generation time: the A>G
signature on the heavy strand that would produce such a shift rises with it across mammals (Mikhailova et al.
2022, NAR 50:10264). Gestation: mtDNA GC content and gestation period were the two most influential traits
in a model of mammalian lifespan (Kulaga et al. 2021, IJMS 22:1073), and in the site's table, gestation and
birth weight follow mammalian GC content on contrasts as closely as lifespan does.

Generation time is not in the build. AnAge's age at female maturity is, and serves as its proxy (`--maturity
any` takes the male age where the female one is missing). The competitors of lifespan are therefore age at
maturity, gestation (or incubation) period and birth weight, each log10. For each competitor T, on the species
with lifespan, adult mass and T (the same species in all its models), per feature, standardised so that a
slope is in standard deviations of the feature, every term reported (mass included):
  lifespan          feature ~ log10 lifespan
  T                 feature ~ log10 T
  lifespan+mass     feature ~ log10 lifespan + log10 adult mass
  T+mass            feature ~ log10 T + log10 adult mass
  lifespan+T+mass   feature ~ log10 lifespan + log10 T + log10 mass
and, on the species with all of them (competitor "all"):
  mass              feature ~ log10 mass
  lifespan          feature ~ log10 lifespan
  all               feature ~ log10 lifespan + maturity + gestation + birth weight + mass
  pace_size         feature ~ slow-fast axis with size
  pace+mass         feature ~ slow-fast axis beyond size + log10 mass
Both axes are the first phylogenetic principal component (Revell 2009, on the correlation matrix), oriented so
that lifespan loads positively and scaled to unit SD across species. "With size" is the component of log10
lifespan, maturity, gestation, birth weight and adult mass. "Beyond size" is the component of the first four
after each is regressed on log10 mass under Brownian motion on the tree.

Body mass is not only a nuisance here: size is part of the slow life history, and large animals live long
partly for the same reasons. A model with mass answers the narrower question of what a trait adds beyond size,
and removes real signal if size lies on the same path; the models without mass give the total association.
Both are reported.

These traits are not independent: a slow life history is long gestation, late maturity, heavy young and long
life together, and an mtDNA feature can be tied to all of them for one evolutionary reason. So the joint
models describe how a signal is shared, and are never read as one trait explaining it "instead of" another
(the trait measured with less error takes the shared part, and maximum lifespan is a noisy record). The
slow-fast axis asks the joint question directly.
Each is a phylogenetic regression on the dated tree (pgls.py), with Pagel's lambda estimated per feature and
model ("pagel") and fixed at 1 ("brownian", the model under which the regression equals the one on
phylogenetically independent contrasts, so it can be compared with the site's contrast test). Fourfold
composition is not a covariate here, unlike H1 and H8: composition is what is being explained.

The features were fixed before any result was seen: GC content, GC skew and AT skew of the protein-coding
genes and of the whole record; A, C, G and T at fourfold-degenerate sites (`fourfold`); the share of
pyrimidine-ending codons that end in C over the protein-coding genes (`c3_share`); and the 32 codons ending
in U or C, each as a percentage of all complete codons, as the site reports them. q: Benjamini-Hochberg over
these features, within class, competitor, tree, model and term.
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
TRAIT_COLUMNS = ["ll", "lm", "lmat", "lgest", "lbw", "pace", "pace_size"]
NAMES = {"ll": "lifespan", "lmat": "maturity", "lgest": "gestation", "lbw": "birth_weight", "lm": "mass", "pace": "pace", "pace_size": "pace_size"}
PACE_TRAITS = ["ll", "lmat", "lgest", "lbw"]
COMPETITORS = ["lmat", "lgest", "lbw"]


def traits(species: pd.DataFrame, maturity: str = "female") -> pd.DataFrame:
    """log10 lifespan, adult mass, age at maturity, gestation (or incubation) and birth weight.
    `maturity` 'any': the male age at maturity where the female one is missing."""
    log = lambda s: np.log10(s.where(s > 0))
    mat = species.female_maturity_days
    if maturity == "any":
        mat = mat.fillna(species.male_maturity_days)
    elif maturity != "female":
        raise ValueError(f"maturity must be 'female' or 'any', not {maturity!r}")
    return pd.DataFrame({"slug": species.slug, "ll": log(species.max_longevity_years),
                         "lm": log(species.adult_weight_g), "lmat": log(mat),
                         "lgest": log(species.gestation_days), "lbw": log(species.birth_weight_g)})


def models() -> list[tuple[str, str, list[str]]]:
    """(competitor, model, columns): every column is a reported term."""
    out = []
    for c in COMPETITORS:
        n = NAMES[c]
        out += [(n, "lifespan", ["ll"]), (n, n, [c]), (n, "lifespan+mass", ["ll", "lm"]), (n, f"{n}+mass", [c, "lm"]),
                (n, f"lifespan+{n}+mass", ["ll", c, "lm"])]
    out += [("all", "mass", ["lm"]), ("all", "lifespan", ["ll"]), ("all", "all", ["ll", *COMPETITORS, "lm"]),
            ("all", "pace_size", ["pace_size"]), ("all", "pace+mass", ["pace", "lm"])]
    return out


def pace_axis(Y: np.ndarray, mass: np.ndarray | None, C: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """The slow-fast axis: first phylogenetic principal component (Revell 2009, correlation matrix) of the
    columns of Y, each first regressed on log10 mass under Brownian motion on C (lambda 1) unless `mass` is None.
    -> (scores scaled to unit standard deviation across species, loadings, share of variance), oriented so
    that the first column (lifespan) loads positively."""
    Ci = np.linalg.inv(C)
    X = np.ones((len(Y), 1)) if mass is None else np.column_stack([np.ones(len(mass)), mass])
    B = np.linalg.solve(X.T @ Ci @ X, X.T @ Ci @ Y)
    R0 = Y - X @ B
    one = np.ones(len(R0))
    a = (one @ Ci @ R0) / (one @ Ci @ one)
    D = R0 - a
    S = D.T @ Ci @ D / (len(D) - 1)
    sd = np.sqrt(np.diag(S))
    w, v = np.linalg.eigh(S / np.outer(sd, sd))
    v1 = v[:, -1] * (1 if v[0, -1] > 0 else -1)
    scores = (D / sd) @ v1
    return scores / scores.std(), v1, float(w[-1] / w.sum())


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
    """The species of one class with a lifespan, adult mass and a tip in the dated tree, with every trait and
    feature (a competitor may be missing: each test drops the species that lack it)."""
    sp = pd.read_parquet(data / "species.parquet")
    sp = sp[sp["class"] == cls]
    steps = [(f"{cls} in the build", len(sp))]
    sp = sp[default_dataset(sp)]
    steps.append(("default dataset (live; lifespan quality high or acceptable; not estimated)", len(sp)))
    tr = traits(sp, maturity).merge(features(data, fourfold), on="slug", how="left")
    tr = tr.dropna(subset=["ll", "lm"])
    steps.append(("and body mass", len(tr)))
    present = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    tr = tr[tr.slug.isin(present)].reset_index(drop=True)
    steps.append(("and a tip in the dated tree", len(tr)))
    for c in COMPETITORS:
        steps.append((f"  of which with {NAMES[c]}" + (f" ({maturity})" if c == "lmat" else ""), int(tr[c].notna().sum())))
    steps.append(("  of which with all three", int(tr[COMPETITORS].notna().all(axis=1).sum())))
    return tr, steps


def run_competitor(sp: pd.DataFrame, tree_file: Path, cls: str, maturity: str, competitor: str,
                   specs: list[tuple[str, list[str]]]) -> tuple[list[dict], dict]:
    """Every model of one competitor, on the species that have all of its columns."""
    cols = sorted({c for _, cs in specs for c in cs}, key=TRAIT_COLUMNS.index)
    sub = sp.dropna(subset=cols).reset_index(drop=True)
    if len(sub) < MIN_SPECIES:
        print(f"{cls}, {competitor}: {len(sub)} species, fewer than {MIN_SPECIES}: not tested")
        return [], {}
    C, order, _ = class_covariance(tree_file, list(sub.slug))
    sub = sub.set_index("slug").loc[order].reset_index()
    raw = [c for c in cols if not c.startswith("pace")]          # the axes are combinations of these
    vif = whitened_vif(C, sub[raw].to_numpy())
    rcols = [c for c in cols if c != "ll"]
    corr = sub[["ll", *rcols]].corr()
    info = {"class": cls, "maturity": maturity, "competitor": competitor, "n": len(sub),
            **{f"vif_{NAMES[c]}": v for c, v in zip(raw, vif)},
            **{f"r_lifespan_{NAMES[c]}": corr.loc["ll", c] for c in rcols}}
    print(f"{cls}, {competitor}: {len(sub)} species; r with lifespan across species "
          + ", ".join(f"{NAMES[c]} {corr.loc['ll', c]:.2f}" for c in rcols)
          + "; VIF on the tree (lambda 1) " + ", ".join(f"{NAMES[c]} {v:.1f}" for c, v in zip(raw, vif)))
    feats = [c for c in sub.columns if c not in ["slug", *TRAIT_COLUMNS]]
    rows = []
    for f in feats:
        x = sub[f].to_numpy(dtype=float)
        ok = np.isfinite(x)
        if ok.sum() < MIN_SPECIES or np.std(x[ok]) == 0:
            continue
        y = (x[ok] - x[ok].mean()) / x[ok].std()
        Cs = C[np.ix_(ok, ok)]
        for tree, lam in (("pagel", None), ("brownian", 1.0)):
            for model, mcols in specs:
                r = fit_terms(Cs, y, sub.loc[ok, mcols].to_numpy(), mcols, lam)
                for col in mcols:
                    b, se, t, p = r[col]
                    rows.append({"class": cls, "maturity": maturity, "competitor": competitor, "tree": tree,
                                 "model": model, "term": NAMES[col], "feature": f, "n": int(ok.sum()),
                                 "lambda": r["lambda"], "slope_sd": b, "se": se, "t": t, "df": r["df"], "p": p})
    return rows, info


def run_class(data: Path, fourfold: Path, tree_file: Path, cls: str, maturity: str) -> tuple[list[dict], list[dict]]:
    sp, steps = select(data, fourfold, tree_file, cls, maturity)
    for label, n in steps:
        print(f"{cls}: {label}: {n}")
    full = sp.dropna(subset=["lm", *PACE_TRAITS]).reset_index(drop=True)
    pace_info = {}
    sp["pace"] = np.nan
    sp["pace_size"] = np.nan
    if len(full) >= MIN_SPECIES:
        C, order, _ = class_covariance(tree_file, list(full.slug))
        full = full.set_index("slug").loc[order].reset_index()
        sp = sp.set_index("slug")
        for col, traits_, mass, label in (("pace", PACE_TRAITS, full.lm.to_numpy(), "beyond size"),
                                          ("pace_size", [*PACE_TRAITS, "lm"], None, "with size")):
            scores, load, share = pace_axis(full[traits_].to_numpy(), mass, C)
            sp.loc[full.slug, col] = scores
            pace_info |= {f"{col}_share": share, **{f"{col}_loading_{NAMES[c]}": l for c, l in zip(traits_, load)}}
            print(f"{cls}: slow-fast axis {label} on {len(full)} species, {100 * share:.0f}% of the variance; loadings "
                  + ", ".join(f"{NAMES[c]} {l:+.2f}" for c, l in zip(traits_, load)))
        sp = sp.reset_index()
    by_comp: dict[str, list[tuple[str, list[str]]]] = {}
    for comp, model, cols in models():
        by_comp.setdefault(comp, []).append((model, cols))
    rows, infos = [], []
    for comp, specs in by_comp.items():
        r, info = run_competitor(sp, tree_file, cls, maturity, comp, specs)
        rows += r
        if info:
            infos.append({**info, **(pace_info if comp == "all" else {})})
    return rows, infos


def main(data: Path, fourfold: Path, tree_file: Path, out: Path, classes: tuple[str, ...] = ("Mammalia", "Aves"),
         maturity: str = "female") -> pd.DataFrame:
    """Writes generation_time.<maturity>.tsv (every fit) and generation_time.<maturity>.classes.tsv (per class
    and competitor: species, the correlation of each trait with lifespan, and the variance inflation)."""
    out.mkdir(parents=True, exist_ok=True)
    rows, infos = [], []
    for cls in classes:
        r, i = run_class(data, fourfold, tree_file, cls, maturity)
        rows += r
        infos += i
    if not rows:
        raise SystemExit("no class had enough species with lifespan, body mass, a competitor and a tip in the tree")
    res = pd.DataFrame(rows)
    res["q"] = res.groupby(["class", "competitor", "tree", "model", "term"]).p.transform(
        lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index))
    res = res.sort_values(["class", "competitor", "tree", "model", "term", "p"]).reset_index(drop=True)
    res.to_csv(out / f"generation_time.{maturity}.tsv", sep="\t", index=False, float_format="%.5g")
    pd.DataFrame(infos).to_csv(out / f"generation_time.{maturity}.classes.tsv", sep="\t", index=False, float_format="%.4g")
    summary = (res[res.tree.eq("pagel")].assign(sig=lambda d: d.q < 0.05)
               .groupby(["class", "competitor", "model", "term"], sort=False)
               .agg(n=("n", "max"), tests=("p", "size"), q_below_005=("sig", "sum"), smallest_p=("p", "min")).reset_index())
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(summary.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res
