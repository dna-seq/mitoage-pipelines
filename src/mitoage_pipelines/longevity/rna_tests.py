"""H37 part A: the composition of each RNA gene against maximum lifespan under the controls H11 and H27 settled, and the
bar of the entry, applied to the composition tests and to the rates (`rna_rates.py`) alike.

The test (`run_class`). Per taxon group (Mammalia, Aves, Actinopterygii, Reptilia: the within-group control), per
gene and per feature, a phylogenetic regression on the build's dated tree (`aa_pgls.fit`, H1's and H11's) of the
standardised feature on log10 lifespan and log10 adult mass (the body-mass control), twice: with Pagel's lambda
estimated (the primary model) and with lambda fixed at 1, the model under which the regression equals the one on
phylogenetically independent contrasts (the PIC control). On H1's species (the default dataset of the manuscript, body
mass, a tip in the tree, all 13 proteins), and again on the species of H27's `high` confidence level. The 13 protein
genes go through the same code: their effects are the comparison of the bar. The features are the composition features
that apply to a stretch with no codons: GC content, GC skew, AT skew (the build's region statistics) and the shares of
A, C, G and T (the analogue of H1's fourfold shares), each in gene sense. The covariance of a class is that of its
species of the set, a gene using the rows of its own species.

The bar (`judge`), as H37 states it, per taxon group and test (a feature, or the rate):
  (1) the gene survives the controls: the lifespan slope has a Holm-adjusted p below 0.05 over the 24 genes, both with
      lambda estimated and with lambda 1, with the same sign. (A rate is a contrast on cherries, the PIC control by
      construction, with the body-mass contrast in the model.)
  (2) its effect, the absolute partial correlation of lifespan, is at least the median of the 13 protein genes' in
      the same group, the same test and the same species set;
  (3) it holds on the high-confidence subset: Holm-adjusted p below 0.05 there too, with the same sign.
Beside that verdict, `correct_all` and `judge_wide` read the same bar against a correction over all 768 cells of a
species set at once, since the entry's correction is within one test and there are 32 of them; the columns of the
wider reading end in `_all`, and `verdict_all` is the verdict under it.
A gene that passes (1) alone is reported as "not counting". One that passes none is a null, said to be "uninformative"
when its alignment has fewer than WEAK_COLUMNS informative columns among the species of the test. A cell with fewer than
MIN_SPECIES species (MIN_CONTRASTS contrasts for a rate) is "not assessed", as is a passing gene for which fewer than
MIN_PROTEIN protein genes give the median.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from ..genes import MT_GENES
from ..seqio import read_fasta
from . import aa_pgls, rates, species_sets
from .rna_genes import RNA_GENES, WEAK_COLUMNS, file_name, gene_features, informative_columns, read_alignment
from .rna_rates import partial_r, rna_alignments

CLASSES = ["Mammalia", "Aves", "Actinopterygii", "Reptilia"]
FEATURES = ["gc_percent", "gc_skew", "at_skew", "A_share", "C_share", "G_share", "T_share"]
MIN_SPECIES = 25
MIN_PROTEIN = 7
ALPHA = 0.05
SETS = ("default", "high")


# --- the species -----------------------------------------------------------------------------------
def test_species(data: Path, composition: Path, tree_file: Path) -> pd.DataFrame:
    """H1's species with the set flags: not retired, lifespan quality high or acceptable and not estimated, body mass, all
    13 proteins and a lifespan (composition.tsv), a tip of the dated tree. Column `high`: H27's confidence level high."""
    if species_sets.active() is not None:
        raise SystemExit(f"{species_sets.ENV} is set: H37 names its species set itself, unset it")
    cols = ["slug", "record_status", "lifespan_quality", "lifespan_record_type", "data_confidence",
            "lifespan_source", "lifespan_sample_size"]
    sp = pd.read_csv(composition, sep="\t").merge(pd.read_parquet(data / "species.parquet", columns=cols), on="slug", how="left")
    sp = sp[aa_pgls.default_dataset(sp)].dropna(subset=["lm"])
    sp["high"] = species_sets.members(sp, "high")
    from Bio import Phylo
    present = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    return sp[sp.slug.isin(present)].reset_index(drop=True)


def protein_features(data: Path) -> pd.DataFrame:
    """The same composition features of the 13 protein genes, from the build's codon alignments (the gene's bases in its
    own sense, gaps removed; the stop codon is not in the alignment)."""
    rows = []
    for g in MT_GENES:
        for slug, seq in read_fasta(data / "alignments" / f"{g}.codon.fasta.gz").items():
            rows.append({"gene": g, "slug": slug, **gene_features(seq)})
    return pd.DataFrame(rows)


# --- the regressions -------------------------------------------------------------------------------
def fit_cell(C: np.ndarray, y: np.ndarray, ll: np.ndarray, lm: np.ndarray, lam: float | None) -> dict:
    """One regression: the standardised feature on log10 lifespan and log10 mass; the slope of lifespan."""
    z = (y - y.mean()) / y.std(ddof=1)
    r = aa_pgls.fit(C, z, np.column_stack([ll, lm]), lam)
    r["r"] = partial_r(r["t"], r["df"])
    return r


def run_class(data: Path, composition: Path, tree_file: Path, genes_dir: Path, cls: str, features: pd.DataFrame | None = None,
              rna: list[str] | None = None) -> pd.DataFrame:
    """Every gene (RNA and protein), feature, species set and lambda model of one taxon group. `features`: a long table
    with `gene`, `slug` and the FEATURES columns; read from `rna_features.tsv` and the codon alignments when None.
    `rna`: the RNA genes to test (all 24 when None); the protein genes are always the comparison."""
    rna = RNA_GENES if rna is None else rna
    sp0 = test_species(data, composition, tree_file)
    if features is None:
        features = pd.concat([pd.read_csv(genes_dir / "rna_features.tsv", sep="\t"), protein_features(data)], ignore_index=True)
    ali = rna_alignments(data, genes_dir, rna)
    rows = []
    for name in SETS:
        sp = sp0[(sp0["class"] == cls) & (sp0.high if name == "high" else True)]
        if len(sp) < MIN_SPECIES:
            continue
        C, order, extra = aa_pgls.class_covariance(tree_file, list(sp.slug))
        sp = sp.set_index("slug").loc[order]
        pos = {s: i for i, s in enumerate(order)}
        print(f"{cls} {name}: {len(sp)} species ({extra} other tips under their common ancestor)", flush=True)
        for gene in [*rna, *MT_GENES]:
            f = features[features.gene == gene]
            f = f[f.slug.isin(pos)]
            f = f.assign(i=f.slug.map(pos)).sort_values("i")
            idx = f.i.to_numpy()
            informative = informative_columns(ali[gene], set(f.slug)) if gene in ali else np.nan
            for feat in FEATURES:
                base = {"set": name, "class": cls, "gene": gene, "group": "rna" if gene in ali else "protein",
                        "feature": feat, "n": len(f), "informative_columns": informative}
                y = f[feat].to_numpy(dtype=float)
                ok = np.isfinite(y)
                if ok.sum() < MIN_SPECIES or np.nanstd(y[ok]) == 0:
                    rows += [{**base, "model": m, "n": int(ok.sum())} for m in ("lambda_ml", "lambda_1")]
                    continue
                i = idx[ok]
                Cg = C[np.ix_(i, i)]
                ll, lm = sp.ll.to_numpy()[i], sp.lm.to_numpy()[i]
                for model, lam in (("lambda_ml", None), ("lambda_1", 1.0)):
                    r = fit_cell(Cg, y[ok], ll, lm, lam)
                    rows.append({**base, "n": int(ok.sum()), "model": model, "lambda": r["lambda"], "slope": r["slope"],
                                 "se": r["se"], "t": r["t"], "df": r["df"], "p": r["p"], "r": r["r"]})
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    return correct(res)


def correct(res: pd.DataFrame) -> pd.DataFrame:
    """Holm and Benjamini-Hochberg over the 24 RNA genes within each set, group, feature and model; the protein genes
    are not corrected (they are the comparison, not the claim)."""
    res = res.copy()
    res["p_holm"] = np.nan
    res["q"] = np.nan
    rna = res[(res.group == "rna") & res.p.notna()]
    for _, g in rna.groupby(["set", "class", "feature", "model"]):
        res.loc[g.index, "p_holm"] = multipletests(g.p, method="holm")[1]
        res.loc[g.index, "q"] = multipletests(g.p, method="fdr_bh")[1]
    return res.sort_values(["set", "class", "feature", "model", "group", "gene"], kind="mergesort").reset_index(drop=True)


def main_class(data: Path, composition: Path, tree_file: Path, genes_dir: Path, out: Path, cls: str) -> pd.DataFrame:
    from .rna_rates import require_control
    require_control(data, out)
    res = run_class(data, composition, tree_file, genes_dir, cls)
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / f"rna_tests.{cls}.tsv", sep="\t", index=False, float_format="%.5g")
    print(f"{cls}: {len(res)} rows")
    return res


# --- the bar ---------------------------------------------------------------------------------------
def long_composition(tests: pd.DataFrame) -> pd.DataFrame:
    """The composition tests of `run_class` in the one shape the bar reads: a row per set, group, class, gene and
    feature, with the primary (lambda estimated) result and the PIC-equivalent (lambda 1) beside it."""
    a = tests[tests.model == "lambda_ml"].set_index(["set", "group", "class", "gene", "feature"])
    b = tests[tests.model == "lambda_1"].set_index(["set", "group", "class", "gene", "feature"])
    d = a[["n", "informative_columns", "slope", "p", "p_holm", "r"]].copy()
    d["slope_pic"] = b.slope
    d["p_pic"] = b.p
    d["p_holm_pic"] = b.p_holm
    d = d.reset_index().rename(columns={"feature": "test"})
    d.insert(0, "kind", "composition")
    return d


def long_rates(rna: pd.DataFrame, protein: pd.DataFrame, informative: pd.DataFrame) -> pd.DataFrame:
    """The rates in the same shape: the mass-adjusted contrast test with saturated pairs left out (H7's headline), RNA
    genes with Holm over the 24, protein genes by dS. `informative`: gene and class -> informative columns."""
    r = rna[(rna.normalisation == "d_mass_adj") & ~rna.includes_saturated].copy()
    p = protein[(protein.normalisation == "dS_mass_adj") & ~protein.includes_saturated].copy()
    out = []
    for d, group in ((r, "rna"), (p, "protein")):
        out.append(d.assign(group=group, test="rate", kind="rate", slope_pic=d.slope, p_pic=d.p, p_holm_pic=d.p_holm))
    d = pd.concat(out, ignore_index=True).merge(informative, on=["gene", "class"], how="left")
    return d[["kind", "set", "group", "class", "gene", "test", "n", "informative_columns", "slope", "p", "p_holm", "r",
              "slope_pic", "p_pic", "p_holm_pic"]]


def correct_all(rows: pd.DataFrame) -> pd.DataFrame:
    """Holm and Benjamini-Hochberg over every RNA cell of a species set at once, beside the entry's own correction.

    The entry's bar corrects over the 24 genes within one test, and there are 32 tests (7 composition features and the
    rate, in each of the 4 taxon groups), so a cell that passes there has been corrected for 23 other genes and for
    none of the other 31 tests. These columns say whether it still passes when the whole family of 768 cells is
    corrected together. Each lambda model is corrected on its own, the bar reading both; the protein genes are not
    corrected, being the comparison and not the claim."""
    rows = rows.copy()
    for src, holm, bh in (("p", "p_holm_all", "q_all"), ("p_pic", "p_holm_pic_all", "q_pic_all")):
        if src not in rows.columns:
            raise ValueError(f"correct_all needs the raw p-values in {src!r}")
        rows[holm] = np.nan
        rows[bh] = np.nan
        cells = rows[(rows.group == "rna") & rows[src].notna()]
        for _, g in cells.groupby("set"):
            rows.loc[g.index, holm] = multipletests(g[src], method="holm")[1]
            rows.loc[g.index, bh] = multipletests(g[src], method="fdr_bh")[1]
    return rows


def judge(rows: pd.DataFrame) -> pd.DataFrame:
    """The bar, over `long_composition` and `long_rates` rows stacked. One row per kind, test, class and RNA gene, on the
    default set, with the high set beside it; see the module's text. The cell's own number of species or contrasts is
    `n`; a missing regression (too few species) is a row with `p` missing and is "not assessed"."""
    out = []
    protein = rows[rows.group == "protein"]
    rna = rows[rows.group == "rna"]
    for (kind, test, cls), _ in rna.groupby(["kind", "test", "class"], sort=True):
        cell = rna[(rna.set == "default") & (rna.kind == kind) & (rna.test == test) & (rna["class"] == cls)]
        prot = protein[(protein.kind == kind) & (protein.test == test) & (protein["class"] == cls) & (protein.set == "default")].dropna(subset=["r"])
        median = float(np.median(prot.r.abs())) if len(prot) >= MIN_PROTEIN else np.nan
        high = rna[(rna.set == "high") & (rna.kind == kind) & (rna.test == test) & (rna["class"] == cls)].set_index("gene")
        have = set(cell.gene)
        for gene in RNA_GENES:
            if gene not in have:
                out.append({"kind": kind, "test": test, "class": cls, "gene": gene, "median_protein_r": median,
                            "n_protein_genes": len(prot), "verdict": "not assessed"})
        for d in cell.sort_values("gene").itertuples(index=False):
            h = high.loc[d.gene] if d.gene in high.index else None
            row = {"kind": kind, "test": test, "class": cls, "gene": d.gene, "n": d.n, "informative_columns": d.informative_columns,
                   "slope": d.slope, "p_holm": d.p_holm, "p_holm_pic": d.p_holm_pic, "r": d.r, "median_protein_r": median,
                   "n_protein_genes": len(prot), "survives_controls": None, "effect_at_least_median": None,
                   "holds_on_high": None,
                   "high_n": np.nan if h is None else h.n, "high_slope": np.nan if h is None else h.slope,
                   "high_p_holm": np.nan if h is None else h.p_holm}
            if pd.isna(d.p):
                row["verdict"] = "not assessed"
            else:
                c1 = bool(d.p_holm < ALPHA and d.p_holm_pic < ALPHA and np.sign(d.slope) == np.sign(d.slope_pic))
                c2 = bool(abs(d.r) >= median) if np.isfinite(median) else None
                c3 = bool(h is not None and pd.notna(h.p_holm) and h.p_holm < ALPHA and np.sign(h.slope) == np.sign(d.slope))
                row.update(survives_controls=c1, effect_at_least_median=c2, holds_on_high=c3)
                weak = pd.isna(d.informative_columns) or d.informative_columns < WEAK_COLUMNS
                if c1 and c2 is None:
                    row["verdict"] = "not assessed (no protein-gene median)"
                elif c1 and c2 and c3:
                    row["verdict"] = "counts"
                elif c1:
                    row["verdict"] = "not counting"
                else:
                    row["verdict"] = "null, uninformative" if weak else "null"
            out.append(row)
    return pd.DataFrame(out)


def judge_wide(rows: pd.DataFrame, v: pd.DataFrame) -> pd.DataFrame:
    """The same bar read against `correct_all`'s wider correction, added to `judge`'s table as its own columns.

    Only conditions (1) and (3) can change: the effect against the protein-gene median is the same number, corrected
    for nothing. A cell that counts under the entry's own correction and not here is one the other 31 tests account
    for, and the entry says so of it in those words."""
    for col in ("p_holm_all", "p_holm_pic_all"):
        if col not in rows.columns:
            raise ValueError("judge_wide reads the columns of correct_all; call it first")
    key = ["kind", "test", "class", "gene"]
    rna = rows[rows.group == "rna"]
    d = rna[rna.set == "default"].set_index(key)
    h = rna[rna.set == "high"].set_index(key)
    out = []
    for row in v.to_dict("records"):
        k = (row["kind"], row["test"], row["class"], row["gene"])
        row.update(p_holm_all=np.nan, p_holm_pic_all=np.nan, high_p_holm_all=np.nan,
                   survives_controls_all=None, holds_on_high_all=None, verdict_all="not assessed")
        if k not in d.index:
            out.append(row)
            continue
        dd, hh = d.loc[k], (h.loc[k] if k in h.index else None)
        row.update(p_holm_all=dd["p_holm_all"], p_holm_pic_all=dd["p_holm_pic_all"],
                   high_p_holm_all=np.nan if hh is None else hh["p_holm_all"])
        if pd.notna(dd["p"]):
            c1 = bool(dd["p_holm_all"] < ALPHA and dd["p_holm_pic_all"] < ALPHA
                      and np.sign(dd["slope"]) == np.sign(dd["slope_pic"]))
            c2 = row["effect_at_least_median"]
            c3 = bool(hh is not None and pd.notna(hh["p_holm_all"]) and hh["p_holm_all"] < ALPHA
                      and np.sign(hh["slope"]) == np.sign(dd["slope"]))
            weak = pd.isna(row["informative_columns"]) or row["informative_columns"] < WEAK_COLUMNS
            row.update(survives_controls_all=c1, holds_on_high_all=c3)
            if c1 and c2 is None:
                row["verdict_all"] = "not assessed (no protein-gene median)"
            elif c1 and c2 and c3:
                row["verdict_all"] = "counts"
            elif c1:
                row["verdict_all"] = "not counting"
            else:
                row["verdict_all"] = "null, uninformative" if weak else "null"
        out.append(row)
    return pd.DataFrame(out)


def summarise(v: pd.DataFrame) -> pd.DataFrame:
    """One row per RNA gene over the cells of `judge`: how many cells ended in each verdict, the smallest Holm-adjusted p of the
    cells that were assessed, and the cell where the gene's effect is largest against its protein-gene median
    (`best_effect_over_median`: |r| divided by the median, so 1 is the median and above 1 is larger)."""
    out = []
    for gene in RNA_GENES:
        g = v[v.gene == gene]
        ok = g[g.verdict != "not assessed"]
        row = {"gene": gene, "cells": len(g), "assessed": len(ok)}
        for name in ("counts", "not counting", "null", "null, uninformative"):
            row["n_" + name.replace(", ", "_").replace(" ", "_")] = int((g.verdict == name).sum())
        row["n_not_assessed"] = int(g.verdict.str.startswith("not assessed").sum())
        row["n_survives_controls"] = int(g.survives_controls.fillna(False).astype(bool).sum())
        row["min_p_holm"] = ok.p_holm.min() if len(ok) else np.nan
        ratio = (ok.r.abs() / ok.median_protein_r).replace([np.inf, -np.inf], np.nan).dropna()
        if len(ratio):
            b = ok.loc[ratio.idxmax()]
            row.update(best_cell=f"{b.kind}:{b.test}:{b['class']}", best_effect_over_median=float(ratio.max()), best_effect_r=b.r)
        out.append(row)
    return pd.DataFrame(out)


def main_verdict(out: Path, genes_dir: Path) -> pd.DataFrame:
    """Reads rna_tests.<Class>.tsv, rna_rates.tsv and protein_rates.tsv from `out`; writes rna_verdict.tsv."""
    tests = pd.concat([pd.read_csv(out / f"rna_tests.{c}.tsv", sep="\t") for c in CLASSES if (out / f"rna_tests.{c}.tsv").is_file()],
                      ignore_index=True)
    comp = long_composition(tests)
    inf = tests[(tests.set == "default") & (tests.group == "rna")][["gene", "class", "informative_columns"]].drop_duplicates()
    rates_rows = long_rates(pd.read_csv(out / "rna_rates.tsv", sep="\t"), pd.read_csv(out / "protein_rates.tsv", sep="\t"), inf)
    rates_rows = rates_rows[rates_rows["class"].isin(CLASSES)]
    rows = correct_all(pd.concat([comp, rates_rows], ignore_index=True))
    v = judge_wide(rows, judge(rows))
    v.to_csv(out / "rna_verdict.tsv", sep="\t", index=False, float_format="%.5g")
    summarise(v).to_csv(out / "rna_gene_summary.tsv", sep="\t", index=False, float_format="%.5g")
    print(v.groupby(["kind", "verdict"]).size().to_string())
    print("under the correction over all cells:")
    print(v.groupby(["kind", "verdict_all"]).size().to_string())
    return v
