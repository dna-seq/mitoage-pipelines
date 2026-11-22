"""H32: the difference vector, the random window sets and their sums, the two statistics against hand-worked values,
and, where `download cybrid` has put the records, the blocks' places on the two hosts."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid, cybrid_blocks as cb, cybrid_dna

GB = paths.CYBRID_GB / "records.gb"
needs_records = pytest.mark.skipif(not GB.is_file(), reason="needs `mitoage download cybrid`")
AL = cybrid_dna.aligner()
HOST = "ACGTTGCAAGGCTTAACCGGTTAACGTCAGTCAAGGTTCC"


class TestDifferenceVector:
    def test_a_mismatch_a_deletion_and_an_insertion_each_land_on_their_place(self):
        assert cb.difference_vector(HOST, HOST, AL).sum() == 0
        v = cb.difference_vector(HOST, HOST[:10] + "T" + HOST[11:], AL)
        assert v.sum() == 1 and v[10] == 1
        v = cb.difference_vector(HOST, HOST[:20] + HOST[21:], AL)
        assert v.sum() == 1 and v[19:23].sum() == 1                  # a deletion in a run may sit on a neighbour
        v = cb.difference_vector(HOST, HOST[:12] + "GGG" + HOST[12:], AL)
        assert v.sum() == 3 and len(v) == len(HOST) and v[9:13].sum() == 3

    def test_an_insertion_before_the_first_place_counts_on_the_first(self):
        v = cb.difference_vector(HOST, "TT" + HOST, AL)
        assert v.sum() == 2 and v[0] == 2 and len(v) == len(HOST)


class TestWindows:
    def test_sets_stay_inside_and_do_not_overlap(self):
        s = cb.random_sets(60, [24, 16, 18], 500, np.random.default_rng(0))
        assert s.shape == (500, 3) and (s >= 0).all() and (s + np.array([24, 16, 18]) <= 60).all()
        for row in s:
            spans = sorted(zip(row, [24, 16, 18]))
            assert all(a + n <= b for (a, n), (b, _) in zip(spans, spans[1:]))
        assert len({tuple(r) for r in s}) > 1

    def test_sums_are_the_sums_over_the_windows(self):
        v = np.arange(10)
        starts = np.array([[0, 5], [2, 8]])
        assert cb.window_sums(v, starts, [2, 2]).tolist() == [0 + 1 + 5 + 6, 2 + 3 + 8 + 9]


class TestStatistics:
    points = pd.DataFrame({"donor": ["a", "a", "b", "b", "c", "c"], "complex": ["I", "IV"] * 3,
                           "status": ["not_reduced", "not_reduced", "reduced", "not_reduced", "reduced", "reduced"],
                           "pct": [np.nan, 100.0, 60.0, np.nan, 40.0, 50.0]})

    def test_the_area_and_the_spearman_against_hand_worked_values(self):
        # scores a 1, b 2, c 3: reduced points hold 2, 3, 3 and the others 1, 1, 2. Of the 9 pairs, 8 have the
        # reduced one higher and one is a tie: (8 + 0.5) / 9
        sc = {"a": np.array([1.0, 3.0]), "b": np.array([2.0, 2.0]), "c": np.array([3.0, 1.0])}
        area, rho = cb.statistics(sc, self.points)
        assert area[0] == pytest.approx(8.5 / 9) and area[1] == pytest.approx(0.5 / 9)
        from scipy.stats import spearmanr
        st = self.points[self.points.pct.notna()]
        want = spearmanr(st.pct, [1.0, 2.0, 3.0, 3.0])[0]
        assert rho[0] == pytest.approx(want) and rho[1] == pytest.approx(-want) and want < 0

    def test_a_disagreeing_point_is_left_out_of_the_area_and_a_constant_score_is_a_half(self):
        p = self.points.copy(); p.loc[0, "status"] = "disagree"
        sc = {"a": np.array([5.0, 1.0]), "b": np.array([5.0, 2.0]), "c": np.array([5.0, 3.0])}
        area, rho = cb.statistics(sc, p)
        assert area[0] == 0.5 and np.isnan(rho[0])
        # without a's first point: reduced hold 2, 3, 3, the others 1 and 2; of 6 pairs 5 are higher and one ties
        assert area[1] == pytest.approx(5.5 / 6)


@pytest.fixture(scope="module")
def features():
    return cybrid_dna.read_features(GB)


@needs_records
class TestOnTheRecords:
    def test_the_mouse_blocks_are_its_records_and_the_human_ones_are_placed_by_alignment(self, features):
        m = features[cybrid.GENOMES["mus-musculus"]]
        got = cb.blocks_of(features, "mus-musculus")
        assert got == {k: (a - m["dloop"][1], b - m["dloop"][1]) for k, (a, b) in m["blocks"].items()}
        h = cb.blocks_of(features, "homo-sapiens")
        seq, start = features[cybrid.GENOMES["homo-sapiens"]]["dloop"]
        n = 16569
        places = {k: ((start + a) % n + 1, (start + b - 1) % n + 1) for k, (a, b) in h.items()}
        assert places == {"CSB1": (210, 234), "CSB2": (299, 314), "CSB3": (346, 363)}
        assert seq[h["CSB3"][0]:h["CSB3"][1]] == "TGCCAAACCCCAAAAACA"

    def test_a_host_against_itself_has_no_difference_in_any_track(self, features):
        for v in cb.tracks(features, "mus-musculus", "mus-musculus", AL).values():
            assert v.sum() == 0

    def test_the_tracks_have_the_hosts_lengths(self, features):
        h = features[cybrid.GENOMES["mus-musculus"]]
        t = cb.tracks(features, "mus-musculus", "mus-spretus", AL)
        assert len(t["control region"]) == len(h["dloop"][0])
        assert len(t["rRNA genes"]) == sum(len(h["rrna"][k][0]) for k in ("12S", "16S"))
        assert 0 < t["control region"].sum() and 0 < t["protein genes"].sum()
