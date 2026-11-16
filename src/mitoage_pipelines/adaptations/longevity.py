"""H26: long-lived for body mass as one more trait in H21's design (docs/HYPOTHESES.md), mammals and birds.

This module defines the trait and calls the code of H21 unchanged (`clades`, `pairs`, `selection`, `coupling`,
`run.reading`). Its tables go under one folder, `out/adapt-longevity/` by default:

  traits.tsv, residuals.tsv              `adapt longevity-traits`
  <Class>/pairs.tsv, pair_tests.tsv      `adapt longevity-pairs`           (R2, both definitions)
  <Class>/selection/long_lived/draw<k>/  `adapt longevity-selection-prep`  (R1, R4), then `adapt relax` and
                                         `adapt absrel` with `--out <folder>/<Class> --trait long_lived`
  Mammalia/coupling.tsv                  `adapt longevity-coupling`        (C, on the set H23 prepared)
  readout.tsv                            `adapt longevity-readout`

The trait, as fixed in H26:

- Species with a residual: the class's species of H1's default dataset (`aa_pgls.default_dataset`) with an adult mass
  above 0, a tip of the dated tree and all 12 heavy-strand codon alignments.
- Residual: ordinary least squares of log10 maximum lifespan on log10 adult mass, one fit per class.
- `long_lived`: residual at or above the class's 75th percentile = 1; below the class's median = 0; between the two
  unknown. A species of the class without a residual is unknown.
- `long_lived_in_order`: the same cuts on the residual of a fit with one intercept per order and a common slope, over
  the orders with at least 5 species with a residual; species of smaller orders are unknown.

Unknown is never written as 0, and `how` says why a state is what it is.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from statsmodels.stats.multitest import multipletests

from ..longevity.aa_pgls import default_dataset
from . import coupling, pairs, selection
from .run import _get, reading
from .traits import COLUMNS

CLASSES = ("Mammalia", "Aves")
TRAIT, TRAIT_IN_ORDER = "long_lived", "long_lived_in_order"
UPPER_Q, LOWER_Q = 0.75, 0.5
MIN_ORDER = 5
RESIDUAL_COLUMNS = ["slug", "class", "order", "log10_lifespan", "log10_mass", "residual", "residual_in_order"]


def live(species: pd.DataFrame) -> pd.DataFrame:
    return species[species.record_status.ne("retired")].reset_index(drop=True)


def residuals(species: pd.DataFrame, usable: set[str], cls: str) -> pd.DataFrame:
    """One row per species of the class with a residual: the class-wide one and, where its order has at least
    MIN_ORDER such species, the one within its order (else NaN)."""
    m = species["class"].eq(cls) & default_dataset(species) & species.adult_weight_g.gt(0) & species.slug.isin(usable)
    d = species.loc[m, ["slug", "class", "order"]].copy()
    d["log10_lifespan"] = np.log10(species.loc[m, "max_longevity_years"].to_numpy(float))
    d["log10_mass"] = np.log10(species.loc[m, "adult_weight_g"].to_numpy(float))
    d = d.sort_values("slug", kind="stable").reset_index(drop=True)
    d["residual"] = np.nan
    d["residual_in_order"] = np.nan
    if len(d) >= 3:
        slope, intercept = np.polyfit(d.log10_mass, d.log10_lifespan, 1)
        d["residual"] = d.log10_lifespan - (intercept + slope * d.log10_mass)
    big = d.order.map(d.order.value_counts()).ge(MIN_ORDER) & d.order.notna() & d.order.ne("")
    k = d[big]
    if k.order.nunique() >= 1 and len(k) >= 3:
        # one intercept per order, a common slope: centre both variables within each order, then one slope
        cm = k.log10_mass - k.groupby("order").log10_mass.transform("mean")
        cl = k.log10_lifespan - k.groupby("order").log10_lifespan.transform("mean")
        slope = float((cm * cl).sum() / (cm * cm).sum()) if float((cm * cm).sum()) > 0 else 0.0
        d.loc[big, "residual_in_order"] = cl - slope * cm
    return d[RESIDUAL_COLUMNS]


def state_of(resid: pd.Series) -> tuple[pd.Series, pd.Series]:
    """1 at or above the 75th percentile of the residuals that exist, 0 below their median, unknown between;
    unknown for a missing residual. Returns (state, how)."""
    r = resid.dropna()
    state = pd.Series(pd.NA, index=resid.index, dtype="Int64")
    how = pd.Series("no residual", index=resid.index, dtype=object)
    if r.empty:
        return state, how
    hi, mid = float(r.quantile(UPPER_Q)), float(r.quantile(LOWER_Q))
    has = resid.notna()
    state[has & resid.ge(hi)] = 1
    state[has & resid.lt(mid)] = 0
    how[has & resid.ge(hi)] = "residual at or above the 75th percentile"
    how[has & resid.lt(mid)] = "residual below the median"
    how[has & resid.ge(mid) & resid.lt(hi)] = "residual between the median and the 75th percentile: intermediate, left out"
    return state, how


def classify(species: pd.DataFrame, usable: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(traits table in `traits.COLUMNS`, residuals table). Every live species of the two classes has a row per trait."""
    sp = live(species)
    parts, res = [], []
    for cls in CLASSES:
        r = residuals(sp, usable, cls)
        res.append(r)
        c = sp[sp["class"].eq(cls)].sort_values("slug", kind="stable").reset_index(drop=True)
        by = r.set_index("slug")
        why = pd.Series("no lifespan record of the default dataset, adult mass, tip or full set of alignments", index=c.index)
        for trait, col in ((TRAIT, "residual"), (TRAIT_IN_ORDER, "residual_in_order")):
            st, how = state_of(c.slug.map(by[col]))
            in_r = c.slug.isin(by.index)
            how = how.where(how.ne("no residual"), why.where(~in_r, f"order with fewer than {MIN_ORDER} species with a residual"))
            parts.append(pd.DataFrame({"trait": trait, "slug": c.slug.to_numpy(), "class": cls, "order": c.order.to_numpy(),
                                       "family": c.family.to_numpy(), "state": pd.array(st.to_numpy(), dtype="Int64"),
                                       "how": how.to_numpy()}))
    t = pd.concat(parts, ignore_index=True)[COLUMNS].sort_values(["trait", "slug"], kind="stable").reset_index(drop=True)
    return t, pd.concat(res, ignore_index=True)


def load(data: Path, tree_file: Path):
    species = live(pd.read_parquet(data / "species.parquet"))
    tips = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    dist = pairs.Distances(data, species)
    return species, tips, dist


def main_traits(data: Path, tree_file: Path, out: Path) -> pd.DataFrame:
    species, tips, dist = load(data, tree_file)
    t, r = classify(species, tips & dist.aligned)
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "traits.tsv", sep="\t", index=False)
    r.to_csv(out / "residuals.tsv", sep="\t", index=False, float_format="%.6g")
    s = t.assign(state=t.state.astype("string").fillna("unknown")).groupby(["trait", "class", "state"]).size().unstack(fill_value=0)
    print(s.to_string())
    return t


def read_traits(out: Path) -> pd.DataFrame:
    f = out / "traits.tsv"
    if not f.is_file():
        raise SystemExit(f"{f} is missing: run `mitoage adapt longevity-traits` first")
    return pd.read_csv(f, sep="\t", dtype={"state": "Int64"})


def run_pairs(data: Path, tree_file: Path, out: Path) -> pd.DataFrame:
    """R2 for both definitions in each class; the size term is the pair's difference in log10 adult mass."""
    species, tips, dist = load(data, tree_file)
    tr = read_traits(out)
    usable = tips & dist.aligned
    mass = np.log10(species.set_index("slug").adult_weight_g.where(lambda w: w > 0))
    all_tests = []
    for cls in CLASSES:
        scope = species.slug[species["class"].eq(cls)]
        tables, tests = [], []
        for trait in (TRAIT, TRAIT_IN_ORDER):
            p = pairs.run(tree_file, tr, trait, usable, dist, scope)
            if p.empty:
                tests.append({"class": cls, "trait": trait, "measure": "d_omega", "with_saturated": False, "n_pairs": 0,
                              "status": "not assessed: no pairs"})
                continue
            tables.append(p)
            for m in ("d_omega", "d_dS", "d_dN"):
                for sat in (False, True):
                    tests.append({"class": cls, "trait": trait, **pairs.test(p, m, sat)})
            tests.append({"class": cls, "trait": trait, **pairs.size_adjusted(p, mass, "d_omega")})
        d = out / cls
        d.mkdir(parents=True, exist_ok=True)
        (pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()).to_csv(d / "pairs.tsv", sep="\t", index=False, float_format="%.6g")
        pd.DataFrame(tests).to_csv(d / "pair_tests.tsv", sep="\t", index=False, float_format="%.6g")
        all_tests += tests
    res = pd.DataFrame(all_tests)
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res


def selection_prep(data: Path, tree_file: Path, out: Path, cls: str) -> pd.DataFrame:
    """RELAX and aBSREL inputs of `long_lived` in one class: <out>/<cls>/selection/long_lived/draw00..draw20."""
    if cls not in CLASSES:
        raise SystemExit(f"--class is one of {', '.join(CLASSES)}; got {cls}")
    species, tips, dist = load(data, tree_file)
    tr = read_traits(out)
    codes = species.set_index("slug").transl_table
    scope = set(species.slug[species["class"].eq(cls)])
    return selection.prepare(data, tree_file, tr, TRAIT, tips & dist.aligned, codes, out / cls / "selection" / TRAIT, scope)


def run_coupling(out: Path, work: Path) -> pd.DataFrame:
    """C1 and C2 of `long_lived` in mammals, on the species set and branch lengths H23 prepared under `work`."""
    tr = read_traits(out)
    tr = tr[tr["class"].eq("Mammalia")]
    if not (work / "genes.txt").is_file():
        rows = [{"readout": r, "trait": TRAIT, "status": f"not assessed: {work} holds no prepared set (H23's `adapt coupling-prep`)"}
                for r in ("C1", "C2")]
    else:
        rows = coupling.run(work, work / "bl", tr, TRAIT)
    res = pd.DataFrame(rows).assign(**{"class": "Mammalia"})
    d = out / "Mammalia"
    d.mkdir(parents=True, exist_ok=True)
    res.to_csv(d / "coupling.tsv", sep="\t", index=False, float_format="%.6g")
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res


def positive_control(adapt_readout: Path) -> bool | None:
    """Whether H23's positive control (bats) passed, read from H21-H25's readout; None when that readout is absent."""
    if not adapt_readout.is_file():
        return None
    r = pd.read_csv(adapt_readout, sep="\t")
    p = r[(r.hypothesis == "H23") & (r.readout == "P")]
    if p.empty:
        return None
    return bool(p.reading.iloc[0] == "passes")


def readout(out: Path, adapt_readout: Path) -> pd.DataFrame:
    """One row per readout of H26, Holm within each class over its primary readouts that were assessed, and the
    reading fixed in H26 (H21's, by `run.reading`)."""
    passed = positive_control(adapt_readout)
    rows = []

    def add(cls, rid, trait, primary, r, estimate, predicted):
        rows.append({"hypothesis": "H26", "class": cls, "readout": rid, "trait": trait, "primary": primary, "kind": "test",
                     "status": r.get("status", "not assessed: not run"), "estimate": r.get(estimate, np.nan),
                     "p": r.get("p", np.nan), "predicted": predicted,
                     "detail": {k: v for k, v in r.items() if k not in ("status", "p")}})

    for cls in CLASSES:
        d = out / cls
        sel = d / "selection" / TRAIT
        add(cls, "R1", TRAIT, True, selection.r1(selection.summarise(sel)) if sel.is_dir() else {}, "K", +1)
        pt = pd.read_csv(d / "pair_tests.tsv", sep="\t") if (d / "pair_tests.tsv").is_file() else pd.DataFrame()
        add(cls, "R2", TRAIT, True, _get(pt, trait=TRAIT, measure="d_omega", with_saturated=False), "median", -1)
        add(cls, "R2", TRAIT_IN_ORDER, False, _get(pt, trait=TRAIT_IN_ORDER, measure="d_omega", with_saturated=False), "median", -1)
        add(cls, "R2 dS", TRAIT, False, _get(pt, trait=TRAIT, measure="d_dS", with_saturated=False), "median", -1)
        if cls == "Mammalia":
            cp = pd.read_csv(d / "coupling.tsv", sep="\t") if (d / "coupling.tsv").is_file() else pd.DataFrame()
            for rid in ("C1", "C2"):
                add(cls, rid, TRAIT, True, _get(cp, readout=rid, trait=TRAIT), "estimate", 0)
    res = pd.DataFrame(rows)
    res["p_holm"] = np.nan
    for _, g in res[res.primary & res.status.eq("assessed") & res.p.notna()].groupby("class"):
        res.loc[g.index, "p_holm"] = multipletests(g.p, method="holm")[1]
    readings = []
    for r in res.itertuples():
        if r.readout in ("C1", "C2") and r.status == "assessed" and passed is None:
            readings.append("not read: H23's positive control (bats) has no result yet")
        else:
            readings.append(reading(r, bool(passed)))
    res["reading"] = readings
    res = res.drop(columns="detail").assign(detail=[str(x) for x in res.detail])
    res.to_csv(out / "readout.tsv", sep="\t", index=False, float_format="%.5g")
    with pd.option_context("display.width", 250, "display.max_colwidth", 70):
        print(res.drop(columns="detail").to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res
