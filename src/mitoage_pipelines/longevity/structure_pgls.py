"""H12: the structural features of the mitochondrial proteins against lifespan, class by class.

Two families of features, each fixed before any result was seen (docs/HYPOTHESES.md, H12):

  contacts    complex IV's mito-nuclear contact mismatch, H9's variable taken into every class. Per species, the
              share of the 298 contact columns of COX1-3 (5Z62, `mitoport burial`) where its residue differs from
              the column's modal residue among the class's tested species, gaps left out; the same over the 298
              burial-matched control columns; and their difference ("contact minus control", the test). H9 asked
              this of all species pooled and of mammals; here each class has its own consensus and its own fit.
  hydropathy  GRAVY (Kyte-Doolittle) of each of the 13 proteins and of all 13 together (`total_pc`), as the build
              publishes it in protein_gene_features.parquet.

The species and the model are H1's (aa_pgls.py): the default dataset with body mass, all 13 proteins,
fourfold counts and a tip of the dated tree; log10 maximum lifespan on the standardised feature (a slope is
log10 years per standard deviation of the feature), with log10 adult mass and the fourfold composition (A4,
C4, T4) in the model ("mass+composition", the primary model), and with mass only ("mass", the model H9 used);
a phylogenetic regression on the class's part of the dated tree, Pagel's lambda estimated by maximum
likelihood for each fit. q: Benjamini-Hochberg within class, family and model.

A class or feature that cannot be fitted is written as a row with status "not assessed" and the reason, never
left out: a test that did not run is not a negative.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..mitoport.iface import mismatch
from ..seqio import read_fasta
from .aa_pgls import class_covariance, fit, select

MIN_SPECIES = 25
CLASSES = ("Mammalia", "Aves", "Actinopterygii", "Reptilia")
MODELS = {"mass+composition": ["lm", "A4", "C4", "T4"], "mass": ["lm"]}
FAMILIES = ("contacts", "hydropathy")
COX = ("COX1", "COX2", "COX3")
GRAVY_REGIONS = ("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
                 "total_pc")
CONTACT_VARIABLES = ("contact minus control", "contact", "control")


def contact_columns(sites: pd.DataFrame) -> dict[str, dict[str, np.ndarray]]:
    """{'contact'|'control': {gene: 0-based alignment columns}} from sites.tsv (aln_col is 1-based)."""
    return {grp: {g: sites[(sites.gene == g) & (sites.group == grp)].aln_col.to_numpy() - 1 for g in COX}
            for grp in ("contact", "control")}


def contact_mismatch(alns: dict[str, dict[str, str]], slugs: list[str], cols: dict[str, dict[str, np.ndarray]]) -> pd.DataFrame:
    """Per species of `slugs`: contact and control mismatch, each the column-weighted mean over COX1-3 of the share
    of filled columns that differ from the modal residue among `slugs`, and their difference."""
    out = {}
    for grp in ("contact", "control"):
        num, den = np.zeros(len(slugs)), 0
        for g in COX:
            share, _ = mismatch(alns[g], slugs, cols[grp][g])
            num += share * len(cols[grp][g]); den += len(cols[grp][g])
        out[grp] = num / den
    res = pd.DataFrame(out, index=pd.Index(slugs, name="slug"))
    res["contact minus control"] = res.contact - res.control
    return res[list(CONTACT_VARIABLES)]


def gravy_table(data: Path) -> pd.DataFrame:
    """GRAVY per species (rows: slug) and region (columns 'gravy:<region>'), retired records left out."""
    pgf = pd.read_parquet(data / "protein_gene_features.parquet", columns=["slug", "region", "gravy"])
    live = set(pd.read_parquet(data / "species.parquet", columns=["slug", "record_status"])
               .query("record_status != 'retired'").slug)
    pgf = pgf[pgf.slug.isin(live) & pgf.region.isin(GRAVY_REGIONS)]
    wide = pgf.pivot_table(index="slug", columns="region", values="gravy", aggfunc="first")
    return wide[[r for r in GRAVY_REGIONS if r in wide.columns]].add_prefix("gravy:")


def fit_features(sp: pd.DataFrame, feats: pd.DataFrame, C: np.ndarray, cls: str, family: str) -> list[dict]:
    """Every feature of `feats` (indexed by slug, aligned to `sp` and `C`) under every model."""
    rows = []
    for f in feats.columns:
        x = feats[f].to_numpy(dtype=float)
        ok = np.isfinite(x)
        base = {"class": cls, "family": family, "feature": f}
        if ok.sum() < MIN_SPECIES or np.std(x[ok]) == 0:
            reason = f"{int(ok.sum())} species with a value" if ok.sum() < MIN_SPECIES else "no variation"
            rows += [{**base, "model": m, "n": int(ok.sum()), "status": f"not assessed: {reason}"} for m in MODELS]
            continue
        z = (x[ok] - x[ok].mean()) / x[ok].std()
        Cs = C[np.ix_(ok, ok)]
        for model, covs in MODELS.items():
            r = fit(Cs, sp.ll.to_numpy()[ok], np.column_stack([z, sp.loc[ok, covs].to_numpy()]))
            h = stats.t.ppf(0.975, r["df"]) * r["se"]
            rows.append({**base, "model": model, "n": int(ok.sum()), "status": "tested", "lambda": r["lambda"],
                         "slope_per_sd": r["slope"], "se": r["se"], "ci95_low": r["slope"] - h,
                         "ci95_high": r["slope"] + h, "t": r["t"], "p": r["p"]})
    return rows


def run_class(data: Path, composition: Path, fourfold: Path, tree_file: Path, cls: str, families: tuple[str, ...],
              sites: pd.DataFrame | None, alns: dict | None) -> list[dict]:
    sp, steps = select(data, composition, fourfold, tree_file, cls)
    for label, n in steps:
        print(f"{cls}: {label}: {n}")
    if len(sp) < MIN_SPECIES:
        print(f"{cls}: fewer than {MIN_SPECIES} species, not assessed")
        return [{"class": cls, "family": fam, "feature": "*", "model": m, "n": len(sp),
                 "status": f"not assessed: {len(sp)} species in the default dataset"} for fam in families for m in MODELS]
    C, order, extra = class_covariance(tree_file, list(sp.slug))
    sp = sp.set_index("slug").loc[order].reset_index()
    print(f"{cls}: {len(sp)} species in the test ({extra} other tips under their common ancestor ignored)")
    rows = []
    if "contacts" in families:
        mm = contact_mismatch(alns, list(sp.slug), contact_columns(sites))
        print(f"{cls}: mean mismatch to the class consensus, contact {mm.contact.mean():.4f}, control {mm.control.mean():.4f}")
        rows += fit_features(sp, mm.loc[sp.slug], C, cls, "contacts")
    if "hydropathy" in families:
        rows += fit_features(sp, gravy_table(data).reindex(sp.slug), C, cls, "hydropathy")
    return rows


def main(data: Path, composition: Path, fourfold: Path, tree_file: Path, work: Path, out: Path,
         classes: tuple[str, ...] = CLASSES, families: tuple[str, ...] = FAMILIES) -> pd.DataFrame:
    """Writes structure_pgls.tsv: class, family, feature, model, n, status, lambda, slope per SD, se, 95% interval,
    t, p, q. `work`: the folder of `mitoport burial` (sites.tsv), needed by the contacts family."""
    unknown = set(families) - set(FAMILIES)
    if unknown:
        raise SystemExit(f"unknown families {sorted(unknown)}; known: {', '.join(FAMILIES)}")
    sites = alns = None
    if "contacts" in families:
        if not (work / "sites.tsv").is_file():
            raise SystemExit(f"{work / 'sites.tsv'} is missing: run `mitoage download pdb 5Z62`, `mitoage mitoport contacts` "
                             "and `mitoage mitoport burial` first, or leave the contacts family out (--family hydropathy)")
        sites = pd.read_csv(work / "sites.tsv", sep="\t")
        alns = {g: read_fasta(data / "alignments" / f"{g}.aa.fasta.gz") for g in COX}
        n = contact_columns(sites)
        print("contact columns:", {g: len(v) for g, v in n["contact"].items()},
              "controls:", {g: len(v) for g, v in n["control"].items()})
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for cls in classes:
        rows += run_class(data, composition, fourfold, tree_file, cls, families, sites, alns)
    res = pd.DataFrame(rows)
    tested = res.status.eq("tested")
    res["q"] = np.nan
    for _, idx in res[tested].groupby(["class", "family", "model"]).groups.items():
        res.loc[idx, "q"] = multipletests(res.loc[idx, "p"], method="fdr_bh")[1]
    order = {c: i for i, c in enumerate(classes)}
    res = res.sort_values(["family", "model", "class", "p"], key=lambda s: s.map(order) if s.name == "class" else s,
                          ascending=[True, False, True, True]).reset_index(drop=True)
    res.to_csv(out / "structure_pgls.tsv", sep="\t", index=False, float_format="%.5g")
    show = res[res.model.eq("mass+composition")]
    summary = (show[show.status.eq("tested")].assign(sig=lambda d: d.q < 0.05)
               .groupby(["family", "class"], sort=False).agg(n=("n", "max"), tests=("p", "size"), q_below_005=("sig", "sum"),
                                                              smallest_p=("p", "min"), smallest_q=("q", "min")).reset_index())
    with pd.option_context("display.width", 220):
        print(summary.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
        na = res[~tested]
        if len(na):
            print("\nnot assessed:\n" + na[["class", "family", "feature", "model", "status"]].to_string(index=False))
    return res
