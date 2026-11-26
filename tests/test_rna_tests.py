"""H37 part A: the composition tests of the RNA genes and the bar of the entry (`longevity/rna_tests.py`).

The regression is checked against the closed-form generalised least squares on the real dated tree and the real species of
the downloaded build (the two rRNA genes of its alignments, which are real RNA genes, against the 13 protein genes). The bar is
checked case by case on tables of exactly the shape the tests write: each case changes one number and names the rule it moves.
"""
import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

from mitoage_pipelines import paths
from mitoage_pipelines.genes import MT_GENES
from mitoage_pipelines.longevity import aa_pgls, rna_genes as rg, rna_tests as rt

COMP = paths.LONGEVITY_OUT / "composition.tsv"
needs_build = pytest.mark.skipif(
    not all(p.is_file() for p in (paths.BUILD / "species.parquet", COMP, paths.DATED_TREE)),
    reason="needs the build, `longevity controls` and `dated-tree`")
RRNA = ["rRNA-12S", "rRNA-16S"]


# --- the regression, on the build -----------------------------------------------------------------
@pytest.fixture(scope="module")
def reptiles():
    features = pd.concat([pd.DataFrame([{"gene": g, "slug": s, **rg.gene_features(seq)}
                                        for s, seq in rg.read_alignment(paths.BUILD / "alignments" / f"{g}.nt.fasta.gz").items()])
                          for g in RRNA] + [rt.protein_features(paths.BUILD)], ignore_index=True)
    res = rt.run_class(paths.BUILD, COMP, paths.DATED_TREE, paths.BUILD, "Reptilia", features, rna=RRNA)
    return features, res


@needs_build
class TestRunClass:
    def test_the_slope_is_the_generalised_least_squares_of_the_standardised_feature_on_lifespan_and_mass(self, reptiles):
        features, res = reptiles
        sp = rt.test_species(paths.BUILD, COMP, paths.DATED_TREE)
        cell = sp[sp["class"] == "Reptilia"]
        C, order, _ = aa_pgls.class_covariance(paths.DATED_TREE, list(cell.slug))
        cell = cell.set_index("slug").loc[order]
        f = features[(features.gene == "COX1") & features.slug.isin(order)].set_index("slug")
        have = [s for s in order if s in f.index]
        i = [order.index(s) for s in have]
        y = f.loc[have, "gc_percent"].to_numpy()
        z = (y - y.mean()) / y.std(ddof=1)
        X = np.column_stack([np.ones(len(have)), cell.ll.to_numpy()[i], cell.lm.to_numpy()[i]])
        Ci = np.linalg.inv(C[np.ix_(i, i)])
        beta = np.linalg.solve(X.T @ Ci @ X, X.T @ Ci @ z)
        row = res[(res.set == "default") & (res.gene == "COX1") & (res.feature == "gc_percent") & (res.model == "lambda_1")].iloc[0]
        assert row.n == len(have) and row.slope == pytest.approx(beta[1], rel=1e-8)
        resid = z - X @ beta
        se = np.sqrt(np.diag(np.linalg.inv(X.T @ Ci @ X)) * (resid @ Ci @ resid / (len(z) - 3)))
        assert row.se == pytest.approx(se[1], rel=1e-8) and row.df == len(z) - 3
        assert row.r == pytest.approx(row.t / np.sqrt(row.t ** 2 + row.df))
        ml = res[(res.set == "default") & (res.gene == "COX1") & (res.feature == "gc_percent") & (res.model == "lambda_ml")].iloc[0]
        assert 0 <= ml["lambda"] <= 1 and row["lambda"] == 1.0

    def test_the_table_has_a_row_per_set_gene_feature_and_model_and_correction_only_over_the_rna_genes(self, reptiles):
        _, res = reptiles
        assert set(res.set) == {"default", "high"} and set(res["class"]) == {"Reptilia"}
        assert len(res) == 2 * (len(RRNA) + len(MT_GENES)) * len(rt.FEATURES) * 2
        assert set(res[res.group == "rna"].gene) == set(RRNA) and set(res[res.group == "protein"].gene) == set(MT_GENES)
        assert res[res.group == "protein"][["p_holm", "q"]].isna().all().all()
        for (s, feat, model), g in res[res.group == "rna"].groupby(["set", "feature", "model"]):
            assert np.allclose(g.p_holm, multipletests(g.p, method="holm")[1]) and np.allclose(g.q, multipletests(g.p, method="fdr_bh")[1])

    def test_the_high_set_is_a_subset_and_the_informative_columns_are_counted_over_the_species_of_the_cell(self, reptiles):
        _, res = reptiles
        a = res[(res.gene == "rRNA-12S") & (res.feature == "gc_percent") & (res.model == "lambda_1")].set_index("set")
        assert 0 < a.loc["high", "n"] < a.loc["default", "n"]
        sp = rt.test_species(paths.BUILD, COMP, paths.DATED_TREE)
        ali = rg.read_alignment(paths.BUILD / "alignments" / "rRNA-12S.nt.fasta.gz")
        for name, cell in (("default", sp[sp["class"] == "Reptilia"]), ("high", sp[(sp["class"] == "Reptilia") & sp.high])):
            assert a.loc[name, "informative_columns"] == rg.informative_columns(ali, set(cell.slug))
        assert a.loc["default", "informative_columns"] != rg.informative_columns(ali)       # not the whole alignment's number

    def test_a_feature_with_no_variation_or_too_few_species_is_a_row_without_a_result(self, reptiles):
        features, _ = reptiles
        f = features.copy()
        f.loc[f.gene == "COX1", "gc_percent"] = 40.0                                   # no variation at all
        cell = rt.test_species(paths.BUILD, COMP, paths.DATED_TREE)
        few = cell[cell["class"] == "Reptilia"].slug.iloc[:rt.MIN_SPECIES - 1]
        f = f[~((f.gene == "COX2") & ~f.slug.isin(few))]                               # a gene seen in too few species of the class
        res = rt.run_class(paths.BUILD, COMP, paths.DATED_TREE, paths.BUILD, "Reptilia", f, rna=[])
        a = res[(res.gene == "COX1") & (res.feature == "gc_percent")]
        assert a.p.isna().all() and a.slope.isna().all() and len(a) == 4                  # two sets, two models, rows kept
        b = res[res.gene == "COX2"]
        assert b.p.isna().all() and b.n.max() == rt.MIN_SPECIES - 1 and len(b) == 2 * len(rt.FEATURES) * 2
        assert res[(res.gene == "COX3") & (res.feature == "gc_percent")].p.notna().all()


# --- the correction ------------------------------------------------------------------------------
def test_holm_is_over_the_rna_genes_within_each_set_class_feature_and_model_and_never_the_protein_genes():
    rng = np.random.default_rng(3)
    rows = []
    for s in ("default", "high"):
        for feat in ("gc_percent", "gc_skew"):
            for model in ("lambda_ml", "lambda_1"):
                for g in rg.RNA_GENES:
                    rows.append({"set": s, "class": "Aves", "gene": g, "group": "rna", "feature": feat, "model": model, "p": rng.uniform(0, 0.2)})
                for g in MT_GENES:
                    rows.append({"set": s, "class": "Aves", "gene": g, "group": "protein", "feature": feat, "model": model, "p": 1e-9})
    rows.append({"set": "default", "class": "Aves", "gene": "tRNA-Ala", "group": "rna", "feature": "gc_skew", "model": "lambda_ml", "p": np.nan})
    res = rt.correct(pd.DataFrame(rows).drop_duplicates(["set", "gene", "feature", "model"], keep="first"))
    for _, g in res[res.group == "rna"].groupby(["set", "class", "feature", "model"]):
        ok = g.p.notna()
        by_hand = np.minimum(1.0, np.maximum.accumulate((len(g[ok]) - np.arange(len(g[ok]))) * np.sort(g[ok].p.to_numpy())))
        assert np.allclose(np.sort(g[ok].p_holm.to_numpy()), by_hand)
        assert (g.p_holm[ok] >= g.p[ok]).all() and (g.p_holm[ok] <= 1).all()
    assert res[res.group == "protein"].p_holm.isna().all()


# --- the bar -------------------------------------------------------------------------------------
MEDIAN = 0.20


def table(kind="composition", test="gc_percent", cls="Mammalia", default=None, high=None, protein=None, extra=None):
    """The rows the bar reads: all 24 RNA genes null on both sets, 13 protein genes whose median |r| is MEDIAN, and the
    changes `default` and `high` ({gene: {column: value}}) make to the RNA genes."""
    base = dict(n=300, informative_columns=40.0, slope=0.1, slope_pic=0.1, p=0.5, p_holm=1.0, p_holm_pic=1.0, r=0.05)
    rows = []
    for s, over in (("default", default or {}), ("high", high or {})):
        for g in rg.RNA_GENES:
            rows.append({"kind": kind, "test": test, "set": s, "group": "rna", "class": cls, "gene": g, **base, **over.get(g, {})})
    r = protein if protein is not None else [MEDIAN + 0.01 * (i - 6) for i in range(13)]
    for s in ("default", "high"):
        for g, v in zip(MT_GENES, r):
            rows.append({"kind": kind, "test": test, "set": s, "group": "protein", "class": cls, "gene": g, **base, "r": v * (50 if s == "high" else 1)})
    return pd.concat([pd.DataFrame(rows), pd.DataFrame(extra or [])], ignore_index=True)


PASS = dict(p=1e-6, p_holm=0.001, p_holm_pic=0.002, r=0.30)
HIGH_OK = dict(p=1e-4, p_holm=0.03, slope=0.1)


def verdicts(**kw):
    return rt.judge(table(**kw)).set_index("gene")


def test_a_gene_that_passes_every_condition_counts():
    v = verdicts(default={"tRNA-Ala": PASS}, high={"tRNA-Ala": HIGH_OK})
    assert v.loc["tRNA-Ala", "verdict"] == "counts"
    assert v.loc["tRNA-Ala", ["survives_controls", "effect_at_least_median", "holds_on_high"]].tolist() == [True, True, True]
    assert v.loc["tRNA-Ala", "median_protein_r"] == pytest.approx(MEDIAN) and v.loc["tRNA-Ala", "n_protein_genes"] == 13
    assert (v.drop(index="tRNA-Ala").verdict == "null").all() and len(v) == 24


def test_passing_the_controls_alone_is_not_counting():
    assert verdicts(default={"tRNA-Ala": PASS})["verdict"]["tRNA-Ala"] == "not counting"             # no high-confidence support
    assert verdicts(default={"tRNA-Ala": {**PASS, "r": 0.19}}, high={"tRNA-Ala": HIGH_OK})["verdict"]["tRNA-Ala"] == "not counting"
    assert verdicts(default={"tRNA-Ala": PASS}, high={"tRNA-Ala": {**HIGH_OK, "p_holm": 0.2}})["verdict"]["tRNA-Ala"] == "not counting"
    assert verdicts(default={"tRNA-Ala": PASS}, high={"tRNA-Ala": {**HIGH_OK, "slope": -0.1}})["verdict"]["tRNA-Ala"] == "not counting"
    assert verdicts(default={"tRNA-Ala": PASS}, high={})["verdict"]["tRNA-Ala"] == "not counting"      # the gene absent from the high set


def test_an_effect_equal_to_the_protein_median_is_large_enough_and_the_sign_must_agree_across_models():
    assert verdicts(default={"tRNA-Ala": {**PASS, "r": MEDIAN}}, high={"tRNA-Ala": HIGH_OK})["verdict"]["tRNA-Ala"] == "counts"
    assert verdicts(default={"tRNA-Ala": {**PASS, "r": -0.30, "slope": -0.1, "slope_pic": -0.1}},
                    high={"tRNA-Ala": {**HIGH_OK, "slope": -0.1}})["verdict"]["tRNA-Ala"] == "counts"      # a negative slope counts by its size
    v = verdicts(default={"tRNA-Ala": {**PASS, "slope_pic": -0.1}}, high={"tRNA-Ala": HIGH_OK})
    assert v["verdict"]["tRNA-Ala"] == "null" and not v["survives_controls"]["tRNA-Ala"]


def test_both_lambda_models_must_pass_the_family_wise_threshold():
    assert verdicts(default={"tRNA-Ala": {**PASS, "p_holm_pic": 0.06}}, high={"tRNA-Ala": HIGH_OK})["verdict"]["tRNA-Ala"] == "null"
    assert verdicts(default={"tRNA-Ala": {**PASS, "p_holm": 0.05}}, high={"tRNA-Ala": HIGH_OK})["verdict"]["tRNA-Ala"] == "null"
    assert verdicts(default={"tRNA-Ala": {**PASS, "p_holm": 0.0499}}, high={"tRNA-Ala": HIGH_OK})["verdict"]["tRNA-Ala"] == "counts"


def test_a_null_on_a_gene_with_few_informative_columns_is_uninformative_not_negative():
    v = verdicts(default={"tRNA-Ala": {"informative_columns": 19.0}, "tRNA-Arg": {"informative_columns": 20.0},
                          "tRNA-Asn": {"informative_columns": float("nan")}})["verdict"]
    assert (v["tRNA-Ala"], v["tRNA-Arg"], v["tRNA-Asn"]) == ("null, uninformative", "null", "null, uninformative")
    weak_pass = verdicts(default={"tRNA-Ala": {**PASS, "informative_columns": 5.0}})
    assert weak_pass["verdict"]["tRNA-Ala"] == "not counting"      # a gene that passes the controls is not a null, however few columns


def test_a_cell_with_no_regression_and_a_gene_with_no_row_are_not_assessed():
    v = verdicts(default={"tRNA-Ala": {"p": np.nan, "p_holm": np.nan, "p_holm_pic": np.nan, "r": np.nan, "slope": np.nan, "n": 12}})
    assert v["verdict"]["tRNA-Ala"] == "not assessed"
    t = table()
    t = t[~((t.set == "default") & (t.gene == "tRNA-Ala"))]
    v = rt.judge(t).set_index("gene")
    assert v["verdict"]["tRNA-Ala"] == "not assessed" and len(v) == 24


def test_the_median_needs_seven_protein_genes_and_comes_from_the_same_group_test_and_set_only():
    few = table(default={"tRNA-Ala": PASS}, high={"tRNA-Ala": HIGH_OK}, protein=[0.1] * 6 + [np.nan] * 7)
    assert rt.judge(few).set_index("gene")["verdict"]["tRNA-Ala"] == "not assessed (no protein-gene median)"
    seven = table(default={"tRNA-Ala": PASS}, high={"tRNA-Ala": HIGH_OK}, protein=[0.1] * 7 + [np.nan] * 6)
    assert rt.judge(seven).set_index("gene")["verdict"]["tRNA-Ala"] == "counts"
    # protein genes of another group, another test or the high set do not move the median of this cell
    noise = pd.concat([table(default={"tRNA-Ala": PASS}, high={"tRNA-Ala": HIGH_OK}),
                       table(cls="Aves", protein=[0.9] * 13), table(test="gc_skew", protein=[0.9] * 13)], ignore_index=True)
    v = rt.judge(noise)
    cell = v[(v["class"] == "Mammalia") & (v.test == "gc_percent")].set_index("gene")
    assert cell.loc["tRNA-Ala", "median_protein_r"] == pytest.approx(MEDIAN) and cell.loc["tRNA-Ala", "verdict"] == "counts"
    assert v.groupby(["kind", "test", "class"]).gene.nunique().eq(24).all() and len(v) == 3 * 24


def test_a_rate_is_judged_by_the_same_bar_and_kept_apart_from_the_composition_tests():
    t = pd.concat([table(kind="composition", default={"tRNA-Ala": PASS}, high={"tRNA-Ala": HIGH_OK}),
                   table(kind="rate", test="rate", default={"tRNA-Ala": {**PASS, "r": 0.1}}, high={"tRNA-Ala": HIGH_OK})], ignore_index=True)
    v = rt.judge(t).set_index(["kind", "gene"])
    assert v.loc[("composition", "tRNA-Ala"), "verdict"] == "counts" and v.loc[("rate", "tRNA-Ala"), "verdict"] == "not counting"


# --- the shapes the bar reads --------------------------------------------------------------------
def test_the_long_tables_put_the_lambda_one_result_beside_the_primary_and_pick_the_headline_rate():
    cols = dict(informative_columns=33.0, n=100, df=97, se=0.1, t=2.0, r=0.2, q=0.5, **{"lambda": 0.5})
    tests = pd.DataFrame([
        {"set": "default", "group": "rna", "class": "Aves", "gene": "tRNA-Ala", "feature": "gc_percent", "model": "lambda_ml", "slope": 0.3, "p": 0.01, "p_holm": 0.2, **cols},
        {"set": "default", "group": "rna", "class": "Aves", "gene": "tRNA-Ala", "feature": "gc_percent", "model": "lambda_1", "slope": -0.2, "p": 0.04, "p_holm": 0.7, **cols}])
    c = rt.long_composition(tests)
    assert len(c) == 1 and c.iloc[0].kind == "composition" and c.iloc[0].test == "gc_percent"
    assert (c.iloc[0].slope, c.iloc[0].slope_pic, c.iloc[0].p_holm, c.iloc[0].p_holm_pic) == (0.3, -0.2, 0.2, 0.7)
    mk = lambda gene, norm, sat, slope: {"set": "default", "class": "Aves", "gene": gene, "normalisation": norm, "includes_saturated": sat,
                                          "n": 50, "slope": slope, "p": 0.01, "p_holm": 0.02, "r": 0.3}
    rna = pd.DataFrame([mk("tRNA-Ala", "d_mass_adj", False, 1.0), mk("tRNA-Ala", "d_raw", False, 2.0), mk("tRNA-Ala", "d_mass_adj", True, 3.0)])
    prot = pd.DataFrame([mk("COX1", "dS_mass_adj", False, 4.0), mk("COX1", "dN_raw", False, 5.0)])
    inf = pd.DataFrame({"gene": ["tRNA-Ala"], "class": ["Aves"], "informative_columns": [41.0]})
    r = rt.long_rates(rna, prot, inf)
    assert sorted(zip(r.group, r.slope)) == [("protein", 4.0), ("rna", 1.0)]
    assert r[r.group == "rna"].informative_columns.iloc[0] == 41.0 and (r.test == "rate").all()


def test_no_composition_test_is_run_without_the_control(tmp_path):
    with pytest.raises(SystemExit, match="rna-control"):
        rt.main_class(tmp_path, tmp_path / "c", tmp_path / "t", tmp_path, tmp_path, "Aves")
