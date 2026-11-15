"""H11: the traits, features and model of `longevity generation-time`, on the downloaded build."""
import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm

from mitoage_pipelines import paths
from mitoage_pipelines.longevity import generation_time as gt

FOURFOLD = paths.LONGEVITY_OUT / "fourfold.tsv"
needs_build = pytest.mark.skipif(
    not (paths.BUILD / "species.parquet").is_file() or not FOURFOLD.is_file(),
    reason="needs `mitoage download build` and `mitoage longevity fourfold`")


@pytest.fixture(scope="module")
def species():
    return pd.read_parquet(paths.BUILD / "species.parquet")


@pytest.fixture(scope="module")
def feats():
    return gt.features(paths.BUILD, FOURFOLD).set_index("slug")


@pytest.fixture(scope="module")
def mammals():
    sp, _ = gt.select(paths.BUILD, FOURFOLD, paths.DATED_TREE, "Mammalia", "female")
    C, order, _ = gt.class_covariance(paths.DATED_TREE, list(sp.slug))
    return sp.set_index("slug").loc[order].reset_index(), C


@needs_build
class TestTraits:
    def test_any_takes_the_male_age_only_where_the_female_one_is_missing(self, species):
        f, a = gt.traits(species, "female"), gt.traits(species, "any")
        has_f = species.female_maturity_days.gt(0).to_numpy()
        assert np.allclose(f.lmat[has_f], a.lmat[has_f])
        only_m = (~species.female_maturity_days.gt(0) & species.male_maturity_days.gt(0)).to_numpy()
        assert only_m.sum() > 0
        assert np.allclose(a.lmat[only_m], np.log10(species.male_maturity_days[only_m]))
        assert f.lmat[only_m].isna().all()

    def test_a_zero_or_missing_value_is_missing_not_minus_infinity(self, species):
        t = gt.traits(species, "any")
        v = t[["ll", "lm", "lmat"]].to_numpy()
        assert not np.isinf(v).any()
        assert np.isnan(v).any()

    def test_an_unknown_maturity_mode_is_refused(self, species):
        with pytest.raises(ValueError):
            gt.traits(species, "median")


@needs_build
class TestFeatures:
    def test_c3_share_is_the_c_ending_share_of_the_pyrimidine_ending_codons(self, feats):
        c = feats[[f"codon:{k}" for k in gt.PYRIMIDINE_ENDING if k[2] == "C"]].sum(axis=1)
        u = feats[[f"codon:{k}" for k in gt.PYRIMIDINE_ENDING if k[2] == "U"]].sum(axis=1)
        ok = feats.c3_share.notna()
        assert ok.sum() > 2000
        assert np.allclose(feats.c3_share[ok], (c / (c + u))[ok])

    def test_codon_percentages_are_of_all_complete_codons(self, feats):
        cu = pd.read_parquet(paths.BUILD / "codon_usage.parquet")
        cu = cu[cu.region.eq("total_pc")].set_index("slug")
        share = cu[gt.PYRIMIDINE_ENDING].sum(axis=1) / cu[gt.CODONS].sum(axis=1) * 100
        got = feats[[f"codon:{k}" for k in gt.PYRIMIDINE_ENDING]].sum(axis=1)
        assert np.allclose(got.loc[share.index], share)
        assert ((got > 0) & (got < 100)).all()

    def test_fourfold_shares_sum_to_one(self, feats):
        ok = feats.A4.notna()
        assert np.allclose(feats.loc[ok, ["A4", "C4", "G4", "T4"]].sum(axis=1), 1)

    def test_gc_content_is_the_build_value(self, feats, species):
        rs = pd.read_parquet(paths.BUILD / "region_stats.parquet")
        rs = rs[rs.region.eq("total_mtDNA") & rs.record_status.ne("retired")].set_index("slug").gc_percent
        assert np.allclose(feats["gc_percent:total_mtDNA"].loc[rs.index], rs)


@needs_build
class TestModel:
    def test_without_the_tree_it_is_ordinary_least_squares(self, mammals):
        sp, C = mammals
        X = sp[["ll", "lmat", "lm"]].to_numpy()
        y = sp["gc_percent:total_mtDNA"].to_numpy()
        r = gt.fit_terms(C, y, X, ["ll", "lmat", "lm"], lam=0.0)
        ols = sm.OLS(y, sm.add_constant(X)).fit()
        for j, n in enumerate(["ll", "lmat", "lm"], start=1):
            assert r[n][0] == pytest.approx(ols.params[j], rel=1e-8)
            assert r[n][3] == pytest.approx(ols.pvalues[j], rel=1e-6)

    def test_without_the_tree_the_inflation_is_the_textbook_one(self, mammals):
        sp, C = mammals
        X = sp[["ll", "lmat", "lm"]].to_numpy()
        vif = gt.whitened_vif(C, X, lam=0.0)
        for j in range(3):
            others = np.delete(X, j, axis=1)
            r2 = sm.OLS(X[:, j], sm.add_constant(others)).fit().rsquared
            assert vif[j] == pytest.approx(1 / (1 - r2), rel=1e-8)

    def test_every_model_of_a_feature_is_fitted_on_the_same_species(self):
        rows, info = gt.run_class(paths.BUILD, FOURFOLD, paths.DATED_TREE, "Mammalia", "female")
        res = pd.DataFrame(rows)
        assert res.groupby("feature").n.nunique().eq(1).all()
        assert info["n"] == res.n.max()
