"""The steps of H21-H25 and their readout (docs/HYPOTHESES.md). Each step writes its tables under one output folder:

  traits.tsv                       `adapt traits`
  pairs.tsv, pair_tests.tsv        `adapt pairs`        (R2 of H21-H23; L3 of H24)
  composition.tsv, lifespan.tsv, l2.tsv, bat_pairs.tsv     `adapt fits`   (R3; H24 P, L1, L2; H25 B1, B2, B3)
  selection/<trait>/draw<k>/       `adapt selection-prep`, `adapt relax`, `adapt absrel`   (R1, R4; H25 B4)
  coupling/<class>/                `adapt coupling-prep`, `mitoport branch-lengths`, `adapt coupling`   (C, P)
  readout.tsv                      `adapt readout`
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from statsmodels.stats.multitest import multipletests

from ..genes import MT_HEAVY_STRAND
from ..longevity.aa_pgls import default_dataset
from ..longevity.generation_time import features
from ..longevity.thermal import traits as thermal_traits
from . import coupling, fits, pairs, selection

# trait -> (hypothesis, role) for the pairs and the selection tests
PAIR_TRAITS = {"diadromy": "H21", "migration": "H21 secondary", "deep_habitat": "H22", "deep_depth": "H22 sensitivity",
               "marine": "H23", "hibernation": "H23"}
SELECTION_TRAITS = {"diadromy": "H21", "deep_habitat": "H22", "marine": "H23", "hibernation": "H23",
                    "bat_hibernation": "H25"}
COMPOSITION_TRAITS = {"diadromy": "Actinopterygii", "migration": "Actinopterygii", "deep_habitat": "Actinopterygii",
                      "deep_depth": "Actinopterygii", "marine": "Mammalia", "hibernation": "Mammalia"}
COUPLING = {"Actinopterygii": ["diadromy", "deep_habitat"], "Mammalia": ["hibernation", "marine", "bats"]}


def load(data: Path, fb: Path, tree_file: Path, traits_file: Path):
    species = pd.read_parquet(data / "species.parquet")
    species = species[species.record_status.ne("retired")].reset_index(drop=True)
    tr = pd.read_csv(traits_file, sep="\t", dtype={"state": "Int64"})
    tips = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    return species, tr, tips


def h24_species(species: pd.DataFrame, tips: set[str], aligned: set[str]) -> pd.Series:
    """H24's species: mammals of H1's default dataset with an adult mass, a tip and all 12 heavy-strand alignments."""
    m = species["class"].eq("Mammalia") & default_dataset(species) & species.adult_weight_g.gt(0)
    return species.slug[m & species.slug.isin(tips & aligned)]


def run_pairs(data: Path, fb: Path, tree_file: Path, traits_file: Path, out: Path) -> pd.DataFrame:
    species, tr, tips = load(data, fb, tree_file, traits_file)
    dist = pairs.Distances(data, species)
    usable = tips & dist.aligned
    size = thermal_traits(species, fb)[0].set_index("slug")["size"]
    tables, tests = [], []
    jobs = [(t, h, usable) for t, h in PAIR_TRAITS.items()]
    jobs.append(("hibernation", "H24 L3", set(h24_species(species, tips, dist.aligned))))
    for trait, hyp, use in jobs:
        p = pairs.run(tree_file, tr, trait, use, dist).assign(hypothesis=hyp)
        tables.append(p)
        measures = ["d_dS", "d_omega", "d_dN"] if hyp == "H24 L3" else ["d_omega", "d_dS", "d_dN"]
        for m in measures:
            for sat in (False, True):
                tests.append({"hypothesis": hyp, "trait": trait, **pairs.test(p, m, sat)})
        tests.append({"hypothesis": hyp, "trait": trait, **pairs.size_adjusted(p, size, measures[0])})
    out.mkdir(parents=True, exist_ok=True)
    pd.concat(tables, ignore_index=True).to_csv(out / "pairs.tsv", sep="\t", index=False, float_format="%.6g")
    res = pd.DataFrame(tests)
    res.to_csv(out / "pair_tests.tsv", sep="\t", index=False, float_format="%.6g")
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res


def run_fits(data: Path, fb: Path, tree_file: Path, traits_file: Path, fourfold: Path, out: Path) -> None:
    species, tr, tips = load(data, fb, tree_file, traits_file)
    feats = features(data, fourfold)
    th = thermal_traits(species, fb)[0].set_index("slug")
    aligned = pairs.Distances(data, species).aligned
    comp = []
    for trait, cls in COMPOSITION_TRAITS.items():
        t = tr[tr.trait.eq(trait) & tr.state.notna() & tr.slug.isin(tips)]
        d = pd.DataFrame({"slug": t.slug, "trait": t.state.astype(float).to_numpy(), "size": t.slug.map(th["size"]).to_numpy()})
        comp.append(fits.feature_fits(d, feats, ["trait", "size"], tree_file, f"R3 {trait}").assign(trait=trait))
    # H24: mammals of the default dataset
    hib = tr[tr.trait.eq("hibernation")].set_index("slug").state
    s24 = h24_species(species, tips, aligned)
    sp = species.set_index("slug")
    d24 = pd.DataFrame({"slug": s24, "ll": np.log10(s24.map(sp.max_longevity_years)), "lm": np.log10(s24.map(sp.adult_weight_g)),
                        "hib": s24.map(hib).astype(float)})
    known24 = d24.dropna(subset=["hib"])
    life = fits.one_fit(known24, "ll", ["hib", "lm"], tree_file, "H24 P")
    comp.append(fits.feature_fits(known24, feats, ["hib", "lm"], tree_file, "H24 L1").assign(trait="hibernation"))
    common = known24.merge(feats, on="slug").slug
    k = known24[known24.slug.isin(common)]
    without = fits.feature_fits(k, feats, ["ll", "lm"], tree_file, "H24 L2 without hibernation")
    with_h = fits.feature_fits(k, feats, ["ll", "lm", "hib"], tree_file, "H24 L2 with hibernation")
    comp += [without.assign(trait="hibernation"), with_h.assign(trait="hibernation")]
    l2 = fits.l2(without, with_h)
    # H25: bats
    bats = species[species.order.eq("Chiroptera") & species.slug.isin(tips & aligned)]
    bh = tr[tr.trait.eq("bat_hibernation")].set_index("slug").state
    d25 = pd.DataFrame({"slug": bats.slug, "ll": np.log10(bats.max_longevity_years).where(default_dataset(bats)),
                        "lm": np.log10(bats.adult_weight_g.where(bats.adult_weight_g > 0)),
                        "hib": bats.slug.map(bh).astype(float)})
    life += fits.one_fit(d25, "ll", ["hib", "lm"], tree_file, "H25 B1")
    comp.append(fits.feature_fits(d25.drop(columns="hib"), feats, ["ll", "lm"], tree_file, "H25 B3 lifespan").assign(trait="bats"))
    comp.append(fits.feature_fits(d25.drop(columns="ll"), feats, ["hib", "lm"], tree_file, "H25 B3 hibernation").assign(trait="bat_hibernation"))
    bp, b2 = fits.bat_pairs(d25.drop(columns="hib").assign(**{"class": "Mammalia"}), pairs.Distances(data, species), tree_file)
    life += b2
    out.mkdir(parents=True, exist_ok=True)
    pd.concat(comp, ignore_index=True).to_csv(out / "composition.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(life).to_csv(out / "lifespan.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame([l2]).to_csv(out / "l2.tsv", sep="\t", index=False, float_format="%.6g")
    bp.to_csv(out / "bat_pairs.tsv", sep="\t", index=False, float_format="%.6g")
    print(pd.DataFrame(life).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(l2)


def selection_prep(data: Path, fb: Path, tree_file: Path, traits_file: Path, out: Path, trait: str) -> pd.DataFrame:
    species, tr, tips = load(data, fb, tree_file, traits_file)
    aligned = pairs.Distances(data, species).aligned
    codes = species.set_index("slug").transl_table
    scope = set(species.slug[species.order.eq("Chiroptera")]) if trait == "bat_hibernation" else None
    return selection.prepare(data, tree_file, tr, trait, tips & aligned, codes, out / "selection" / trait, scope)


def coupling_prep(data: Path, nucaln: Path, tree_file: Path, traits_file: Path, out: Path, cls: str) -> dict:
    species = pd.read_parquet(data / "species.parquet")
    tr = pd.read_csv(traits_file, sep="\t", dtype={"state": "Int64"})
    trait = COUPLING[cls][0]
    ts = set(tr.slug[tr.trait.eq(trait) & tr.state.eq(1)])
    return coupling.prepare(data, nucaln, tree_file, species, cls, ts, out / "coupling" / cls)


def run_coupling(traits_file: Path, out: Path, cls: str) -> pd.DataFrame:
    tr = pd.read_csv(traits_file, sep="\t", dtype={"state": "Int64"})
    work = out / "coupling" / cls
    rows = []
    for trait in COUPLING[cls]:
        rows += coupling.run(work, work / "bl", tr, trait, positive_control=trait == "bats")
    res = pd.DataFrame(rows).assign(**{"class": cls})
    res.to_csv(work / "coupling.tsv", sep="\t", index=False, float_format="%.6g")
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res


# --- the readout -------------------------------------------------------------------------------------------------
def _get(df: pd.DataFrame, **kw) -> dict:
    d = df
    for k, v in kw.items():
        if k not in d:
            return {}
        d = d[d[k].eq(v)]
    return d.iloc[0].to_dict() if len(d) else {}


def readout(out: Path) -> pd.DataFrame:
    """One row per pre-registered readout of H21-H25, with Holm within each hypothesis and the reading fixed there."""
    pt = pd.read_csv(out / "pair_tests.tsv", sep="\t") if (out / "pair_tests.tsv").is_file() else pd.DataFrame()
    life = pd.read_csv(out / "lifespan.tsv", sep="\t") if (out / "lifespan.tsv").is_file() else pd.DataFrame()
    l2 = pd.read_csv(out / "l2.tsv", sep="\t") if (out / "l2.tsv").is_file() else pd.DataFrame()
    coup = pd.concat([pd.read_csv(f, sep="\t") for f in (out / "coupling").glob("*/coupling.tsv")], ignore_index=True) \
        if (out / "coupling").is_dir() and list((out / "coupling").glob("*/coupling.tsv")) else pd.DataFrame()
    rows = []
    def add(h, rid, trait, primary, r, estimate=None, predicted=0, kind="test"):
        st = r.get("status", "not assessed: not run")
        rows.append({"hypothesis": h, "readout": rid, "trait": trait, "primary": primary, "kind": kind, "status": st,
                     "estimate": r.get(estimate) if estimate else np.nan, "p": r.get("p", np.nan),
                     "predicted": predicted, "detail": {k: v for k, v in r.items() if k not in ("status", "p")}})
    for trait, h, pred in (("diadromy", "H21", +1), ("deep_habitat", "H22", 0), ("marine", "H23", 0), ("hibernation", "H23", 0),
                           ("bat_hibernation", "H25", 0)):
        d = out / "selection" / trait
        r1 = selection.r1(selection.summarise(d)) if d.is_dir() else {}
        add(h, "B4" if h == "H25" else "R1", trait, True, r1, "K", pred)
    for trait, h, pred in (("diadromy", "H21", -1), ("deep_habitat", "H22", 0), ("marine", "H23", 0), ("hibernation", "H23", 0)):
        add(h, "R2", trait, True, _get(pt, hypothesis=h, trait=trait, measure="d_omega", with_saturated=False), "median", pred)
    add("H24", "L3", "hibernation", True, _get(pt, hypothesis="H24 L3", measure="d_dS", with_saturated=False), "median", -1)
    for trait, h in (("diadromy", "H21"), ("hibernation", "H23"), ("marine", "H23")):
        for rid in ("C1", "C2"):
            add(h, rid, trait, True, _get(coup, readout=rid, trait=trait), "estimate", 0)
    for rid in ("C1", "C2"):     # H22: outside its correction, computed if the fish set holds enough deep-sea species
        add("H22", rid, "deep_habitat", False, _get(coup, readout=rid, trait="deep_habitat"), "estimate", 0)
    add("H23", "P", "bats", False, _get(coup, readout="P", trait="bats"), "estimate", +1, kind="positive control")
    add("H24", "P", "hibernation", False, _get(life, model="H24 P", term="hib"), "slope", +1, kind="positive control")
    add("H25", "B1", "bat_hibernation", False, _get(life, model="H25 B1", term="hib"), "slope", +1, kind="positive control")
    add("H25", "B2", "bats", True, _get(life, model="B2", response="log10_dS_per_myr", term="ll"), "slope", 0)
    l2r = l2.iloc[0].to_dict() if len(l2) else {}
    rows.append({"hypothesis": "H24", "readout": "L2", "trait": "hibernation", "primary": True, "kind": "reading",
                 "status": l2r.get("status", "not assessed: not run"), "estimate": l2r.get("median_ratio", np.nan),
                 "p": np.nan, "predicted": 0, "detail": l2r})
    res = pd.DataFrame(rows)
    res["p_holm"] = np.nan
    for h, g in res[res.primary & res.kind.eq("test") & res.status.eq("assessed") & res.p.notna()].groupby("hypothesis"):
        res.loc[g.index, "p_holm"] = multipletests(g.p, method="holm")[1]
    p_ok = res[(res.readout == "P") & (res.hypothesis == "H23")]
    p_pass = bool(len(p_ok) and p_ok.status.iloc[0] == "assessed" and p_ok.p.iloc[0] < 0.05 and p_ok.estimate.iloc[0] > 0)
    res["reading"] = [reading(r, p_pass) for r in res.itertuples()]
    res = res.drop(columns="detail").assign(detail=[str(d) for d in res.detail])
    res.to_csv(out / "readout.tsv", sep="\t", index=False, float_format="%.5g")
    with pd.option_context("display.width", 250, "display.max_colwidth", 70):
        print(res.drop(columns="detail").to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res


def reading(r, positive_control_passed: bool) -> str:
    if r.status != "assessed":
        return r.status
    if r.kind == "reading":
        return getattr(r, "detail", {}).get("reading", "") if isinstance(r.detail, dict) else ""
    if r.kind == "positive control":
        return "passes" if (r.p < 0.05 and r.estimate > 0) else "fails"
    if r.readout in ("C1", "C2") and not positive_control_passed:
        return "not read: H23's positive control (bats) did not pass"
    p = r.p_holm if r.primary and not np.isnan(r.p_holm) else r.p
    if r.readout in ("R1", "B4"):
        ctl = r.detail.get("controls_below", 0) if isinstance(r.detail, dict) else 0
        if p >= 0.05:
            return "not detected"
        if ctl < 19:
            return "not separated from the control clades"
        word = "intensified" if r.estimate > 1 else "relaxed"
        if r.predicted:
            return f"{word}: " + ("as predicted" if np.sign(r.estimate - 1) == r.predicted else "opposite to the prediction")
        return word
    if p >= 0.05:
        return "not detected"
    sign = "positive" if r.estimate > 0 else "negative"
    if r.predicted:
        return f"{sign}: " + ("as predicted" if np.sign(r.estimate) == r.predicted else "opposite to the prediction")
    return sign
