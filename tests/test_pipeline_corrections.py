"""Regression cases for missing data and tied observations in the analysis pipelines."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from mitoage_pipelines.adaptations import coupling, pairs
from mitoage_pipelines.longevity import rates, structure_pgls
from mitoage_pipelines.mitoport import cybrid


@pytest.mark.parametrize("missing", ["---", "GCN", "TAA"])
def test_sisters_use_the_same_codons_against_the_outgroup(missing):
    a, b, o = "GCT" * 40 + "GCC" * 20, "GCT" * 40 + missing * 20, "GCT" * 60
    # Pairwise deletion invents a difference although the sisters agree everywhere B is known.
    assert rates.distance(a, o, 2)["dS"] > 0
    first, second = rates.outgroup_distances(a, b, o, 2)
    assert first == second == {"dS": 0.0, "dN": 0.0, "codons": 40}


def test_the_shared_codon_floor_applies_to_both_sisters():
    a, b, o = "GCT" * 40, "GCT" * 20 + "---" * 20, "GCT" * 40
    first, second = rates.outgroup_distances(a, b, o, 2)
    assert first["codons"] == second["codons"] == 20
    assert np.isnan(first["dS"]) and np.isnan(second["dS"])


def test_h7_builds_contrasts_on_the_shared_mask(monkeypatch, tmp_path):
    aln = {"a": "GCT" * 40 + "GCC" * 20, "b": "GCT" * 40 + "---" * 20, "o": "GCT" * 60}
    monkeypatch.setattr(rates, "MT_GENES", ("X",))
    monkeypatch.setattr(rates, "read_fasta", lambda _: aln)
    tree = tmp_path / "tree.nwk"
    tree.write_text("((a:1,b:1):1,o:2);\n")
    sp = pd.DataFrame({"slug": ["a", "b", "o"], "class": "Mammalia", "transl_table": 2,
                       "ll": [1., 2., 3.], "lm": [1., 2., 3.]})
    contrasts = rates.build_contrasts(tmp_path, sp, tree)
    assert set(contrasts.gene) == {"X", rates.CONCAT}
    assert contrasts.codons.eq(40).all() and contrasts.d_dS.eq(0).all()


def test_adaptation_distances_cache_the_shared_mask_in_sister_order(monkeypatch, tmp_path):
    aln = {"a": "GCT" * 38 + "GCC" * 22, "b": "GCT" * 40 + "---" * 20, "o": "GCT" * 60}
    monkeypatch.setattr(pairs, "read_fasta", lambda _: aln)
    sp = pd.DataFrame({"slug": ["a", "b", "o"], "transl_table": 2})
    dist = pairs.Distances(tmp_path, sp, genes=("X",))
    first, second = dist.to_outgroup("a", "b", "o")
    assert first["codons"] == second["codons"] == 40
    assert first["dS"] > second["dS"] == 0
    assert dist.to_outgroup("b", "a", "o") == (second, first)
    assert len(dist.triplet_cache) == 1
    row = pairs.pair_rates([{"taken": True, "sp_trait": "a", "sp_other": "b", "outgroup": "o"}], dist).iloc[0]
    assert row.dS_trait == first["dS"] and row.dS_other == second["dS"]


def test_rate_ties_leave_the_sign_test_but_stay_in_the_regression():
    d = pd.DataFrame({"d_ll": 1., "d_dS": [-1.] * 5 + [1.] * 5 + [0.] * 10})
    result = rates.test_contrasts(d, "d_dS")
    assert result["n"] == 20 and result["sign_n"] == 10
    assert result["lower_rate_in_longer_lived"] == 5
    assert result["sign_p"] == stats.binomtest(5, 10).pvalue == 1.


def test_an_all_tied_rate_has_no_sign_test():
    result = rates.test_contrasts(pd.DataFrame({"d_ll": [1.] * 20, "d_dS": [0.] * 20}), "d_dS")
    assert result["sign_n"] == 0 and np.isnan(result["sign_p"])


@pytest.mark.parametrize("missing", [np.nan, np.inf, -np.inf])
def test_c2_missing_observed_contrast_is_not_assessed(missing):
    result = coupling.c2(pd.Series([missing, 1., 2.]), np.array([True, False, False]),
                         [np.array([False, True, False])] * 1000)
    assert result["status"].startswith("not assessed") and "p" not in result


def test_c2_uses_only_finite_null_draws_and_requires_a_hundred():
    logm = pd.Series([np.nan, 1., 2., 3.])
    fg = np.array([False, True, False, False])
    bad = np.array([True, False, False, False])
    good = np.array([False, False, True, False])
    result = coupling.c2(logm, fg, [bad] + [good] * 100)
    assert result["status"] == "assessed" and result["null_draws"] == 100
    assert result["invalid_null_draws"] == 1 and np.isfinite(result["p"])
    assert coupling.c2(logm, fg, [bad] * 1000 + [good] * 99)["status"].startswith("not assessed")


def test_c2_a_real_zero_is_assessed():
    fg = np.array([True, False, False])
    result = coupling.c2(pd.Series([1., 1., 1.]), fg, [fg] * 100)
    assert result["status"] == "assessed" and result["estimate"] == 0 and result["p"] == 1


def scores():
    return pd.DataFrame([{"host": "host", "donor": donor, "complex": cx, "score": value}
                         for donor, value in (("good", 1.), ("bad", 2.)) for cx in cybrid.STRUCTURES])


@pytest.mark.parametrize("missing", [np.nan, np.inf, -np.inf])
def test_cybrid_total_requires_every_finite_complex_score(missing):
    data = scores()
    data.loc[data.donor.eq("good") & data.complex.eq("IV"), "score"] = missing
    totals = cybrid.complete_totals(data, ("score",)).set_index("donor").score
    assert np.isnan(totals["good"]) and totals["bad"] == 8


def test_cybrid_total_requires_the_complex_row_to_exist():
    data = scores().query("donor != 'good' or complex != 'IV'")
    totals = cybrid.complete_totals(data, ("score",)).set_index("donor").score
    assert np.isnan(totals["good"])
    # H18 legitimately has three modelled complexes for the mouse, supplied explicitly.
    data = scores().query("complex != 'V'")
    totals = cybrid.complete_totals(data, ("score",), {"host": ("I", "III", "IV")}).set_index("donor").score
    assert totals.to_dict() == {"good": 3., "bad": 6.}


@pytest.mark.parametrize("missing_donor", ["good", "bad"])
def test_t4_does_not_read_incomplete_totals_as_support(missing_donor):
    data = scores()
    data = data[~(data.donor.eq(missing_donor) & data.complex.eq("IV"))]
    measurements = pd.DataFrame([
        {"host": "host", "donor": "good", "complex": "I", "outcome": "measured", "significant": "no",
         "pct": 100., "source_doi": "10.test"},
        {"host": "host", "donor": "bad", "complex": "all", "outcome": "no_respiring_cybrid",
         "significant": "not_stated", "pct": np.nan, "source_doi": "10.test"},
    ])
    _, summary = cybrid.summarise(data, measurements, names=("score",), t2b=(), t3=(), t4=("score",))
    result = summary[summary.test.eq("T4")].iloc[0]
    assert np.isnan(result.value) and "not assessed" in result.reads


def test_structure_main_writes_an_entirely_unassessed_table(monkeypatch, tmp_path):
    monkeypatch.setattr(structure_pgls, "select", lambda *args: (pd.DataFrame({"slug": range(10)}), []))
    result = structure_pgls.main(Path("unused"), Path("unused"), Path("unused"), Path("unused"),
                                 tmp_path, tmp_path, classes=("Mammalia",), families=("hydropathy",))
    assert len(result) == len(structure_pgls.MODELS)
    assert result.status.str.startswith("not assessed").all()
    assert result.p.isna().all() and result.q.isna().all()
    saved = pd.read_csv(tmp_path / "structure_pgls.tsv", sep="\t")
    assert len(saved) == len(result) and saved.status.str.startswith("not assessed").all()
