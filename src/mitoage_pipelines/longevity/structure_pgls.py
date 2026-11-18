"""H12: the structural features of the mitochondrial proteins against lifespan, class by class.

Two families of features, each fixed before any result was seen (docs/HYPOTHESES.md, H12):

  contacts    complex IV's mito-nuclear contact mismatch, H9's variable taken into every class. Per species, the
              share of the 298 contact columns of COX1-3 (5Z62, `mitoport burial`) where its residue differs from
              the column's modal residue among the class's tested species, gaps left out; the same over the 298
              burial-matched control columns; and their difference ("contact minus control", the test). H9 asked
              this of all species pooled and of mammals; here each class has its own consensus and its own fit.
  hydropathy  GRAVY (Kyte-Doolittle) of each of the 13 proteins and of all 13 together (`total_pc`), as the build
              publishes it in protein_gene_features.parquet.
              Also fitted under a third model, "mass+own composition": the shares of A, C and U at the fourfold
              sites of the same region (gene sense, from the same table) in place of the 12-gene composition, which
              leaves out ND6, the one protein read from the other strand.
  embedding   (only when asked for, after `structure embed` on a GPU) H9's question in ESM-2 space: per species and
              gene, the cosine distance of its mean ESM-2 650M embedding over the contact residues from the class's
              mean vector, the same over the control residues, weighted over COX1-3 by column count, and their
              difference (the test). A species needs all three genes.

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
OPTIONAL_FAMILIES = ("embedding",)          # needs `structure embed`; run only when named
OWN = "mass+own composition"
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


def own_composition(data: Path) -> dict[str, pd.DataFrame]:
    """{region: shares a, c, u of A, C and U at its fourfold sites, gene sense, per slug}, retired records left out."""
    pgf = pd.read_parquet(data / "protein_gene_features.parquet", columns=["slug", "region", "ff_sites", "ff_a", "ff_c", "ff_u"])
    live = set(pd.read_parquet(data / "species.parquet", columns=["slug", "record_status"])
               .query("record_status != 'retired'").slug)
    pgf = pgf[pgf.slug.isin(live) & pgf.region.isin(GRAVY_REGIONS)]
    out = {}
    for r, d in pgf.groupby("region"):
        n = d.ff_sites.where(d.ff_sites > 0)
        out[r] = pd.DataFrame({"a": d.ff_a / n, "c": d.ff_c / n, "u": d.ff_u / n}).set_axis(d.slug)
    return out


def fit_own(sp: pd.DataFrame, feats: pd.DataFrame, own: dict[str, pd.DataFrame], C: np.ndarray, cls: str) -> list[dict]:
    """Every GRAVY feature under the model with its own region's fourfold composition (and mass)."""
    rows = []
    for f in feats.columns:
        cov = own[f.split(":", 1)[1]].reindex(sp.slug)[["a", "c", "u"]].to_numpy(dtype=float)
        x = feats[f].to_numpy(dtype=float)
        ok = np.isfinite(x) & np.isfinite(cov).all(axis=1)
        base = {"class": cls, "family": "hydropathy", "feature": f, "model": OWN, "n": int(ok.sum())}
        if ok.sum() < MIN_SPECIES or np.std(x[ok]) == 0:
            rows.append({**base, "status": f"not assessed: {int(ok.sum())} species with a value and its composition"})
            continue
        z = (x[ok] - x[ok].mean()) / x[ok].std()
        r = fit(C[np.ix_(ok, ok)], sp.ll.to_numpy()[ok], np.column_stack([z, sp.lm.to_numpy()[ok], cov[ok]]))
        h = stats.t.ppf(0.975, r["df"]) * r["se"]
        rows.append({**base, "status": "tested", "lambda": r["lambda"], "slope_per_sd": r["slope"], "se": r["se"],
                     "ci95_low": r["slope"] - h, "ci95_high": r["slope"] + h, "t": r["t"], "p": r["p"]})
    return rows


def embedding_divergence(emb: pd.DataFrame, slugs: list[str], weights: dict[str, int]) -> pd.DataFrame:
    """Per species of `slugs` with all three genes: cosine distance of its mean embedding over the contact (control)
    residues from the mean of those vectors over `slugs`, weighted over COX1-3 by `weights`; and the difference."""
    dims = [c for c in emb.columns if c.startswith("e") and c[1:].isdigit()]
    emb = emb[emb.slug.isin(slugs)]
    out = {}
    for grp in ("contact", "control"):
        num, have = pd.Series(0.0, index=pd.Index(slugs, name="slug")), pd.Series(0, index=pd.Index(slugs, name="slug"))
        for g in COX:
            d = emb[(emb.gene == g) & (emb.set == grp)].set_index("slug")
            X = d[dims].to_numpy(dtype=np.float64)
            c = X.mean(axis=0)
            dist = 1 - X @ c / (np.linalg.norm(X, axis=1) * np.linalg.norm(c))
            num = num.add(pd.Series(dist * weights[g], index=d.index), fill_value=np.nan)
            have = have.add(pd.Series(1, index=d.index), fill_value=0)
        out[grp] = (num / sum(weights[g] for g in COX)).where(have == len(COX))
    res = pd.DataFrame(out)
    res["contact minus control"] = res.contact - res.control
    return res[list(CONTACT_VARIABLES)]


def run_class(data: Path, composition: Path, fourfold: Path, tree_file: Path, cls: str, families: tuple[str, ...],
              sites: pd.DataFrame | None, alns: dict | None, emb: pd.DataFrame | None = None) -> list[dict]:
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
        g = gravy_table(data).reindex(sp.slug)
        rows += fit_features(sp, g, C, cls, "hydropathy")
        rows += fit_own(sp, g, own_composition(data), C, cls)
    if "embedding" in families:
        n = contact_columns(sites)["contact"]
        dv = embedding_divergence(emb, list(sp.slug), {g: len(v) for g, v in n.items()})
        print(f"{cls}: {int(dv.contact.notna().sum())} species with embeddings of COX1-3; mean distance to the class mean, "
              f"contact {dv.contact.mean():.4f}, control {dv.control.mean():.4f}")
        rows += fit_features(sp, dv.reindex(sp.slug), C, cls, "embedding")
    return rows


def main(data: Path, composition: Path, fourfold: Path, tree_file: Path, work: Path, out: Path,
         classes: tuple[str, ...] = CLASSES, families: tuple[str, ...] = FAMILIES, embed: Path | None = None) -> pd.DataFrame:
    """Writes structure_pgls.tsv: class, family, feature, model, n, status, lambda, slope per SD, se, 95% interval,
    t, p, q. `work`: the folder of `mitoport burial` (sites.tsv), needed by the contacts and embedding families;
    `embed`: the table of `structure embed`, needed by the embedding family."""
    unknown = set(families) - set(FAMILIES) - set(OPTIONAL_FAMILIES)
    if unknown:
        raise SystemExit(f"unknown families {sorted(unknown)}; known: {', '.join(FAMILIES + OPTIONAL_FAMILIES)}")
    sites = alns = emb = None
    if "embedding" in families:
        if embed is None or not embed.is_file():
            raise SystemExit(f"the embedding family needs the table of `mitoage structure embed` ({embed}), made on a GPU")
        emb = pd.read_parquet(embed)
        if not (work / "sites.tsv").exists():
            raise SystemExit(f"{work / 'sites.tsv'} is missing: run `mitoage mitoport burial` first")
        sites = pd.read_csv(work / "sites.tsv", sep="\t")
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
        rows += run_class(data, composition, fourfold, tree_file, cls, families, sites, alns, emb)
    res = pd.DataFrame(rows).reindex(columns=["class", "family", "feature", "model", "n", "status", "lambda",
                                            "slope_per_sd", "se", "ci95_low", "ci95_high", "t", "p"])
    tested = res.status.eq("tested")
    res["q"] = np.nan
    for _, idx in res[tested].groupby(["class", "family", "model"]).groups.items():
        res.loc[idx, "q"] = multipletests(res.loc[idx, "p"], method="fdr_bh")[1]
    order = {c: i for i, c in enumerate(classes)}
    res = res.sort_values(["family", "model", "class", "p"], key=lambda s: s.map(order) if s.name == "class" else s,
                          ascending=[True, False, True, True]).reset_index(drop=True)
    res.to_csv(out / "structure_pgls.tsv", sep="\t", index=False, float_format="%.5g")
    show = res[res.model.isin(["mass+composition", OWN])]
    summary = (show[show.status.eq("tested")].assign(sig=lambda d: d.q < 0.05)
               .groupby(["family", "model", "class"], sort=False).agg(n=("n", "max"), tests=("p", "size"), q_below_005=("sig", "sum"),
                                                              smallest_p=("p", "min"), smallest_q=("q", "min")).reset_index())
    with pd.option_context("display.width", 220):
        print(summary.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
        na = res[res.status.ne("tested")]
        if len(na):
            print("\nnot assessed:\n" + na[["class", "family", "feature", "model", "status"]].to_string(index=False))
    return res
