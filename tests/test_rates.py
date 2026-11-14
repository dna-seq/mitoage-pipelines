"""H7: the counting, the cherries and the contrast test of `longevity rates`."""
import io

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo

from mitoage_pipelines.longevity import rates


def tree_of(newick: str):
    return Phylo.read(io.StringIO(newick), "newick")


class TestGeneticCode:
    def test_fourfold_codons_have_one_synonymous_site_at_the_third_position(self):
        _, syn = rates.code_tables(2)
        for codon in ("GCT", "GCC", "GCA", "GCG"):          # alanine, a fourfold family
            assert syn[codon] == pytest.approx(1.0)

    def test_a_twofold_codon_has_a_third_of_a_site_per_synonymous_change(self):
        aa, syn = rates.code_tables(2)
        # TTT/TTC are Phe, TTA/TTG are Leu in the vertebrate code: one of three changes is synonymous
        assert aa["TTT"] == "F" and aa["TTA"] == "L"
        assert syn["TTT"] == pytest.approx(1 / 3)

    def test_stop_codons_carry_no_synonymous_sites_and_are_not_counted(self):
        aa, syn = rates.code_tables(2)
        stops = [c for c, a in aa.items() if a == "*"]
        assert stops and all(syn[c] == 0 for c in stops)
        assert list(rates.counted_codons(stops[0], stops[0], 2)) == []

    def test_the_vertebrate_and_invertebrate_codes_differ_where_they_should(self):
        vert, _ = rates.code_tables(2)
        invert, _ = rates.code_tables(5)
        assert vert["AGA"] == "*" and invert["AGA"] == "S"


class TestCounting:
    def test_identical_sequences_have_no_differences_and_zero_distance(self):
        s = "GCTGCCGCAGCG"
        assert rates.differences(s, s, 2) == (0.0, 0.0)
        S, N, codons = rates.sites(s, s, 2)
        assert codons == 4
        assert S + N == pytest.approx(3 * codons)

    def test_a_third_position_change_in_a_fourfold_codon_is_one_synonymous_difference(self):
        Sd, Nd = rates.differences("GCT", "GCC", 2)
        assert (Sd, Nd) == (1.0, 0.0)

    def test_a_first_position_change_that_alters_the_amino_acid_is_nonsynonymous(self):
        aa, _ = rates.code_tables(2)
        assert aa["GCT"] != aa["CCT"]
        assert rates.differences("GCT", "CCT", 2) == (0.0, 1.0)

    def test_a_codon_pair_differing_twice_is_averaged_over_its_pathways(self):
        Sd, Nd = rates.differences("GCT", "CCC", 2)       # Ala -> Pro, two changes
        assert Sd + Nd == pytest.approx(2.0)              # every path makes exactly two changes
        assert 0 < Sd < 2 and 0 < Nd < 2                  # neither path is wholly one kind

    def test_a_pathway_through_a_stop_codon_is_not_taken(self):
        aa, _ = rates.code_tables(2)
        # TCA (Ser) to TAC (Tyr) differs at two positions. One order passes through TAA, a stop in the
        # vertebrate code, and must be refused; the other passes through TCC (Ser).
        assert aa["TAA"] == "*" and aa["TCC"] == "S" and aa["TCA"] == "S" and aa["TAC"] == "Y"
        Sd, Nd = rates.differences("TCA", "TAC", 2)
        assert (Sd, Nd) == (1.0, 1.0)          # the single legal path: synonymous then nonsynonymous
        # averaging both paths, as a version that ignored stops would, gives 0.5 and 1.5 instead
        assert (Sd, Nd) != (0.5, 1.5)

    def test_gaps_and_ambiguous_codons_are_skipped_entirely(self):
        assert list(rates.counted_codons("GCT---GCC", "GCTGCCGCC", 2)) == [("GCT", "GCT"), ("GCC", "GCC")]
        assert list(rates.counted_codons("GCN", "GCT", 2)) == []

    def test_distance_refuses_to_report_below_the_codon_floor(self):
        short = "GCT" * (rates.MIN_CODONS - 1)
        assert np.isnan(rates.distance(short, short, 2)["dS"])
        enough = "GCT" * rates.MIN_CODONS
        assert rates.distance(enough, enough, 2)["dS"] == 0.0


class TestJukesCantor:
    def test_zero_divergence_is_zero_and_saturation_is_unknown_not_a_number(self):
        assert rates.jukes_cantor(0.0) == 0.0
        assert np.isnan(rates.jukes_cantor(0.75))
        assert np.isnan(rates.jukes_cantor(0.9))

    def test_the_correction_raises_the_observed_proportion_and_is_monotone(self):
        ps = [0.05, 0.1, 0.2, 0.4, 0.6]
        ds = [rates.jukes_cantor(p) for p in ps]
        assert all(d > p for d, p in zip(ds, ps))
        assert ds == sorted(ds)


class TestCherries:
    def test_only_pairs_of_usable_tips_that_share_a_parent_are_cherries(self):
        t = tree_of("(((a:1,b:1):1,c:2):1,(d:1,e:1):2);")
        found = {tuple(sorted((x, y))) for _, x, y in rates.cherries(t, {"a", "b", "c", "d", "e"})}
        assert found == {("a", "b"), ("d", "e")}

    def test_a_cherry_whose_sister_is_unusable_is_not_offered(self):
        t = tree_of("(((a:1,b:1):1,c:2):1,(d:1,e:1):2);")
        found = {tuple(sorted((x, y))) for _, x, y in rates.cherries(t, {"a", "c", "d", "e"})}
        assert found == {("d", "e")}

    def test_a_polytomy_is_not_a_cherry(self):
        t = tree_of("((a:1,b:1,c:1):1,d:2);")
        assert list(rates.cherries(t, {"a", "b", "c", "d"})) == []

    def test_the_outgroup_is_the_nearest_usable_tip_outside_the_cherry(self):
        # two tips sit outside the cherry under the same ancestor: near is closer than mid
        t = tree_of("(((a:1,b:1):1,(near:0.5,mid:9):0.5):5,far:20);")
        parents = rates.parent_map(t)
        node = next(n for n, x, y in rates.cherries(t, {"a", "b"}))
        usable = {"a", "b", "near", "mid", "far"}
        assert rates.outgroup_for(t, parents, node, {"a", "b"}, usable) == "near"
        assert t.distance(node, "near") < t.distance(node, "mid")

    def test_the_search_climbs_when_nothing_nearer_is_usable(self):
        t = tree_of("(((a:1,b:1):1,c:0.5):5,far:20);")
        parents = rates.parent_map(t)
        node = next(n for n, x, y in rates.cherries(t, {"a", "b"}))
        assert rates.outgroup_for(t, parents, node, {"a", "b"}, {"a", "b", "far"}) == "far"
        assert rates.outgroup_for(t, parents, node, {"a", "b"}, {"a", "b"}) is None


class TestContrastTest:
    def make(self, n=40, slope=0.5, seed=0):
        rng = np.random.default_rng(seed)
        d_ll = rng.normal(size=n)
        return pd.DataFrame({"d_ll": d_ll, "d_lm": rng.normal(size=n), "time_my": rng.uniform(1, 10, n),
                             "d_dS": slope * d_ll, "d_dN": rng.normal(size=n)})

    def test_a_planted_slope_is_recovered_through_the_origin(self):
        r = rates.test_contrasts(self.make(slope=0.5), "d_dS")
        assert r["slope"] == pytest.approx(0.5, abs=1e-9)
        assert r["n"] == 40

    def test_a_fit_with_no_residual_reports_no_t_and_no_p_rather_than_certainty(self):
        # the planted data above is noiseless: the slope is exact, so its standard error is zero
        r = rates.test_contrasts(self.make(slope=0.5), "d_dS")
        assert r["se"] == pytest.approx(0.0, abs=1e-12)
        assert np.isnan(r["t"]) and np.isnan(r["p"])

    def test_a_noisy_slope_does_get_a_t_and_a_p(self):
        df = self.make(slope=0.5)
        df["d_dS"] = df.d_dS + np.random.default_rng(1).normal(scale=0.2, size=len(df))
        r = rates.test_contrasts(df, "d_dS")
        assert np.isfinite(r["t"]) and 0 < r["p"] < 1

    def test_no_intercept_is_fitted(self):
        # contrasts carry no intercept: a constant offset must not be absorbed
        df = self.make(slope=0.0)
        df["d_dS"] = 1.0
        assert rates.test_contrasts(df, "d_dS")["slope"] != pytest.approx(0.0, abs=1e-6)

    def test_the_sign_test_counts_pairs_where_the_longer_lived_sister_is_slower(self):
        df = pd.DataFrame({"d_ll": [1.0] * 10 + [-1.0] * 10, "d_dS": [-1.0] * 10 + [1.0] * 10,
                           "time_my": 1.0, "d_lm": 0.0})
        r = rates.test_contrasts(df, "d_dS")
        assert r["lower_rate_in_longer_lived"] == 20 and r["sign_n"] == 20
        flipped = df.assign(d_dS=-df.d_dS)
        assert rates.test_contrasts(flipped, "d_dS")["lower_rate_in_longer_lived"] == 0

    def test_pairs_with_no_lifespan_difference_leave_the_sign_test(self):
        df = pd.DataFrame({"d_ll": [1.0, -1.0, 0.0] * 10, "d_dS": [-1.0, 1.0, 1.0] * 10,
                           "time_my": 1.0, "d_lm": 0.0})
        r = rates.test_contrasts(df, "d_dS")
        assert r["sign_n"] == 20 and r["n"] == 30

    def test_dividing_by_time_changes_the_slope_it_should(self):
        df = self.make(slope=0.5)
        df["time_my"] = 2.0
        raw = rates.test_contrasts(df, "d_dS")["slope"]
        per = rates.test_contrasts(df, "d_dS", per_time=True)["slope"]
        assert per == pytest.approx(raw / 2)

    def test_too_few_contrasts_is_no_answer_rather_than_a_weak_one(self):
        assert rates.test_contrasts(self.make(n=rates.MIN_CONTRASTS - 1), "d_dS") is None
        assert rates.test_contrasts(self.make(n=rates.MIN_CONTRASTS), "d_dS") is not None

    def test_rows_without_the_values_the_model_needs_are_dropped(self):
        df = self.make(n=30)
        df.loc[:4, "d_lm"] = np.nan
        assert rates.test_contrasts(df, "d_dS")["n"] == 30
        assert rates.test_contrasts(df, "d_dS", covariates=["d_lm"])["n"] == 25
