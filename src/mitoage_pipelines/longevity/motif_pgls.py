"""H8: the mtDNA motif, repeat and contact-zone features of the database against lifespan.

The features are the ones the site's pipeline writes into `mtdna_structure.parquet` - G-quadruplex
motifs, homopurine and homopyrimidine tracts, homopolymer runs, short tandem repeats, direct,
inverted and mirror repeats, and the repeat pairs bridging the two deletion contact zones. They are
read from the build, never recomputed here, so what is tested is what the database publishes.

The model is H1's, unchanged (`pgls.py`): log10 maximum lifespan on the feature, with log10 adult
body mass and the fourfold-site composition, as a phylogenetic regression on the build's dated tree
with Pagel's lambda estimated per class. Features are standardised, so a slope is the change in
log10 lifespan per standard deviation of the feature. q is Benjamini-Hochberg within class and
subset.

Two subsets are reported where the data allow it: every record, and records whose sequencing
technology is known to include long reads. mitoVGP and MitoHiFi both report that short-read
reference mitogenomes misrepresent exactly these features, so the second subset is the one that
would settle whether a correlation is biology or assembly. RefSeq mitochondrion records almost never
state their technology, so that subset is usually empty and the test says so rather than quietly
reporting only the first.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import optimize, stats
from statsmodels.stats.multitest import multipletests

from .aa_pgls import class_covariance, default_dataset
from .pgls import gls, nll

MIN_SPECIES = 25
COVARIATES = ["lm", "A4", "C4", "T4"]
# columns of mtdna_structure.parquet that are coordinates or keys rather than features
NOT_FEATURES = {"species_id", "slug", "refseq", "orientation", "repeat_min_length",
                "major_arc_start", "major_arc_end", "contact_zone1_start", "contact_zone1_end",
                "contact_zone2_start", "contact_zone2_end"}
LONG_READ = ("pacbio", "nanopore", "ont", "hifi", "minion", "promethion", "gridion", "sequel", "revio")


def feature_columns(structure: pd.DataFrame) -> list[str]:
    """Numeric columns of the structure table that are features, in the file's own order."""
    return [c for c in structure.columns
            if c not in NOT_FEATURES and pd.api.types.is_numeric_dtype(structure[c])]


def long_read_slugs(structure: pd.DataFrame, technology: Path | None = None) -> set[str]:
    """Species whose record names a long-read technology: from technology.tsv (`longevity technology`)
    when given, else from the build's own column; empty when neither exists."""
    if technology is not None:
        t = pd.read_csv(technology, sep="\t")
        return set(t.slug[t.technology_class == "long_read"])
    if "sequencing_technology" not in structure.columns:
        return set()
    tech = structure.sequencing_technology.fillna("").str.lower()
    return set(structure.slug[tech.apply(lambda t: any(k in t for k in LONG_READ))])


def select(data: Path, composition: Path, fourfold: Path, tree_file: Path, cls: str) -> pd.DataFrame:
    sp = pd.read_csv(composition, sep="\t")
    sp = sp[sp["class"] == cls].merge(
        pd.read_parquet(data / "species.parquet",
                        columns=["slug", "record_status", "lifespan_quality", "lifespan_record_type"]),
        on="slug", how="left")
    sp = sp[default_dataset(sp)].dropna(subset=["lm"])
    sp = sp.merge(pd.read_csv(fourfold, sep="\t"), on="slug")
    present = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    return sp[sp.slug.isin(present)].reset_index(drop=True)


def test_features(sp: pd.DataFrame, structure: pd.DataFrame, features: list[str], tree_file: Path,
                  cls: str, subset: str) -> list[dict]:
    sp = sp.merge(structure, on="slug")
    if len(sp) < MIN_SPECIES:
        return []
    C, order, _ = class_covariance(tree_file, list(sp.slug))
    sp = sp.set_index("slug").loc[order].reset_index()
    if len(sp) < MIN_SPECIES:
        return []
    y = sp.ll.to_numpy()
    base = np.column_stack([np.ones(len(sp)), sp[COVARIATES].to_numpy()])
    lam = float(optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(base, y, C), method="bounded").x)
    rows = []
    for f in features:
        x = pd.to_numeric(sp[f], errors="coerce").to_numpy(dtype=float)
        ok = np.isfinite(x)
        if ok.sum() < MIN_SPECIES or np.std(x[ok]) == 0:
            continue
        z = (x[ok] - x[ok].mean()) / x[ok].std()
        X = np.column_stack([base[ok], z])
        b, se, df = gls(X, y[ok], C[np.ix_(ok, ok)], lam)
        t = b[-1] / se[-1]
        rho, rho_p = stats.spearmanr(x[ok], y[ok])
        rows.append({"class": cls, "subset": subset, "feature": f, "n": int(ok.sum()), "lambda": lam,
                     "spearman_rho": rho, "spearman_p": rho_p,
                     "slope_per_sd": b[-1], "se": se[-1], "t": t, "p": 2 * stats.t.sf(abs(t), df)})
    return rows


def main(data: Path, composition: Path, fourfold: Path, tree_file: Path, out: Path,
         structure_path: Path | None = None, technology: Path | None = None,
         classes: tuple[str, ...] = ("Mammalia", "Aves", "Actinopterygii", "Reptilia")) -> pd.DataFrame:
    """Writes motif_pgls.tsv: class, subset, feature, n, lambda, Spearman, PGLS slope per SD, p, q.

    The species, lifespans and covariates come from `data` (the current build, as every other test
    here reads it); the features come from `structure_path`, which defaults to the same build. They
    are separable because the motif features are not yet part of a released build: until that work is
    merged, the table is pointed at the branch that writes it, and the species the two share are the
    ones tested. The count is printed, never assumed.
    """
    out.mkdir(parents=True, exist_ok=True)
    path = structure_path or (data / "mtdna_structure.parquet")
    if not path.is_file():
        raise SystemExit(
            f"{path} is not in this build. The motif features ship from the site pipeline's "
            "mtdna-features work; point --build at a build that contains them.")
    structure = pd.read_parquet(path)
    features = feature_columns(structure)
    if technology is None and (out / "technology.tsv").is_file():
        technology = out / "technology.tsv"
    long_read = long_read_slugs(structure, technology)
    print(f"sequencing technology read from: {technology or 'the build (no technology.tsv)'}")
    print(f"{path}: {len(structure)} species, {len(features)} numeric features")
    if structure_path is not None:
        build_slugs = set(pd.read_parquet(data / "species.parquet", columns=["slug"]).slug)
        shared = len(build_slugs & set(structure.slug))
        print(f"features from a separate build: {shared} of the current build's {len(build_slugs)} species "
              f"carry them; the rest cannot be tested")
    print(f"records naming a long-read technology: {len(long_read)}"
          + ("" if long_read else "  (the long-read-only subset cannot be tested on this build)"))
    rows = []
    for cls in classes:
        sp = select(data, composition, fourfold, tree_file, cls)
        rows += test_features(sp, structure, features, tree_file, cls, "all_records")
        if long_read:
            rows += test_features(sp[sp.slug.isin(long_read)], structure, features, tree_file, cls, "long_read")
    if not rows:
        raise SystemExit("no class had enough species with both a lifespan and the motif features")
    res = pd.DataFrame(rows)
    res["q"] = res.groupby(["class", "subset"]).p.transform(
        lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index))
    res = res.sort_values("p").reset_index(drop=True)
    res.to_csv(out / "motif_pgls.tsv", sep="\t", index=False, float_format="%.5g")
    with pd.option_context("display.width", 220):
        print(res.head(20).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"\n{len(res)} tests; q < 0.05 in {(res.q < 0.05).sum()}")
    return res
