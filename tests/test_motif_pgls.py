"""H8: the feature selection, the long-read subset and the model of `longevity motif-pgls`."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines.longevity import motif_pgls


def structure(n=60, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "species_id": np.arange(n), "slug": [f"sp-{i}" for i in range(n)],
        "refseq": [f"NC_{i:06d}.1" for i in range(n)], "orientation": "canonical",
        "repeat_min_length": 10, "major_arc_start": 1, "major_arc_end": 100,
        "contact_zone1_start": 1, "contact_zone1_end": 10,
        "contact_zone2_start": 20, "contact_zone2_end": 30,
        "g4_forward_total": rng.integers(0, 5, n), "mirror_repeats_total": rng.integers(100, 900, n),
        "contact_zone_direct_pairs": rng.integers(0, 40, n),
    })


class TestFeatureSelection:
    def test_coordinates_and_keys_are_not_treated_as_features(self):
        cols = motif_pgls.feature_columns(structure())
        assert set(cols) == {"g4_forward_total", "mirror_repeats_total", "contact_zone_direct_pairs"}

    def test_every_coordinate_column_of_the_published_file_is_excluded(self):
        # the zone and arc coordinates describe where a feature was measured, not how much of it there is
        for c in ("major_arc_start", "major_arc_end", "contact_zone1_start", "contact_zone2_end",
                  "repeat_min_length", "species_id"):
            assert c in motif_pgls.NOT_FEATURES
            assert c not in motif_pgls.feature_columns(structure())

    def test_text_columns_are_not_features(self):
        s = structure()
        s["sequencing_technology"] = "Illumina"
        assert "sequencing_technology" not in motif_pgls.feature_columns(s)

    def test_the_file_order_is_kept_so_the_output_is_stable(self):
        s = structure()
        assert motif_pgls.feature_columns(s) == [c for c in s.columns if c in motif_pgls.feature_columns(s)]


class TestLongReadSubset:
    def test_a_build_without_the_column_yields_an_empty_subset_not_an_error(self):
        assert motif_pgls.long_read_slugs(structure()) == set()

    def test_long_read_platforms_are_recognised_whatever_the_case_or_wording(self):
        s = structure(n=6)
        s["sequencing_technology"] = ["PacBio Sequel II", "Oxford Nanopore MinION", "Illumina HiSeq",
                                      "illumina novaseq; pacbio", None, "Sanger dideoxy"]
        assert motif_pgls.long_read_slugs(s) == {"sp-0", "sp-1", "sp-3"}

    def test_short_read_only_and_missing_technologies_stay_out(self):
        s = structure(n=3)
        s["sequencing_technology"] = ["Illumina MiSeq", None, ""]
        assert motif_pgls.long_read_slugs(s) == set()


class TestModel:
    def make(self, n=60, effect=0.0, seed=1):
        rng = np.random.default_rng(seed)
        x = rng.normal(size=n)
        sp = pd.DataFrame({
            "slug": [f"sp-{i}" for i in range(n)],
            "ll": effect * x + rng.normal(scale=0.3, size=n),
            "lm": rng.normal(size=n), "A4": rng.normal(size=n),
            "C4": rng.normal(size=n), "T4": rng.normal(size=n),
        })
        st = pd.DataFrame({"slug": sp.slug, "feature_x": x})
        return sp, st

    def tree_file(self, tmp_path, n=60):
        # a star tree: every tip at the same depth, so the covariance is the identity scaled by that depth
        newick = "(" + ",".join(f"sp-{i}:1" for i in range(n)) + ");"
        p = tmp_path / "star.nwk"
        p.write_text(newick)
        return p

    def test_a_planted_effect_is_found_and_a_null_one_is_not(self, tmp_path):
        tf = self.tree_file(tmp_path)
        sp, st = self.make(effect=1.0)
        strong = motif_pgls.test_features(sp, st, ["feature_x"], tf, "Testaceae", "all_records")
        assert len(strong) == 1 and strong[0]["p"] < 1e-6
        sp0, st0 = self.make(effect=0.0, seed=7)
        weak = motif_pgls.test_features(sp0, st0, ["feature_x"], tf, "Testaceae", "all_records")
        assert weak[0]["p"] > 0.05

    def test_the_slope_is_per_standard_deviation_so_rescaling_the_feature_changes_nothing(self, tmp_path):
        tf = self.tree_file(tmp_path)
        sp, st = self.make(effect=0.8)
        a = motif_pgls.test_features(sp, st, ["feature_x"], tf, "T", "all_records")[0]
        b = motif_pgls.test_features(sp, st.assign(feature_x=st.feature_x * 1000 + 7),
                                     ["feature_x"], tf, "T", "all_records")[0]
        assert a["slope_per_sd"] == pytest.approx(b["slope_per_sd"], rel=1e-9)
        assert a["p"] == pytest.approx(b["p"], rel=1e-9)

    def test_a_feature_present_in_too_few_species_is_no_answer(self, tmp_path):
        # plenty of species, but the feature itself is measured in fewer than the floor
        tf = self.tree_file(tmp_path)
        sp, st = self.make(effect=1.0)
        st.loc[motif_pgls.MIN_SPECIES - 1:, "feature_x"] = np.nan
        assert motif_pgls.test_features(sp, st, ["feature_x"], tf, "T", "all_records") == []

    def test_too_few_species_is_no_answer(self, tmp_path):
        n = motif_pgls.MIN_SPECIES - 1
        tf = self.tree_file(tmp_path, n)
        sp, st = self.make(n=n, effect=1.0)
        assert motif_pgls.test_features(sp, st, ["feature_x"], tf, "T", "all_records") == []

    def test_a_constant_feature_is_skipped_rather_than_divided_by_zero(self, tmp_path):
        tf = self.tree_file(tmp_path)
        sp, st = self.make()
        assert motif_pgls.test_features(sp, st.assign(feature_x=5.0), ["feature_x"], tf, "T", "all_records") == []

    def test_species_missing_the_feature_are_dropped_and_counted(self, tmp_path):
        tf = self.tree_file(tmp_path)
        sp, st = self.make(effect=0.5)
        st.loc[:9, "feature_x"] = np.nan
        r = motif_pgls.test_features(sp, st, ["feature_x"], tf, "T", "all_records")[0]
        assert r["n"] == 50
