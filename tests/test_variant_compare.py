"""H34, part A: joining MitImpact's scores to H29's variants and the table of areas.

Expected values are worked out by hand in each test from the small tables it builds."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines.mitoport import cybrid, variant_compare as vc


def mitimpact_csv(tmp_path, rows):
    cols = ["Start", "Ref", "Alt", "Gene_symbol", "functional_effect_general", "AlphaMissense_score", "ESM1b_score", "APOGEE2_score", "PhyloP_100V", "MLC_score"]
    f = tmp_path / "m.csv"
    pd.DataFrame(rows, columns=cols).to_csv(f, index=False)
    return f


ROWS = [(10, "A", "C", "MT-ND1", "missense_variant", 0.9, -12.0, 0.8, 5.0, 0.7),
        (10, "A", "G", "MT-ND1", "missense_variant", 0.1, -2.0, 0.1, 5.0, 0.7),
        (20, "T", "C", "MT-ATP8", "missense_variant", 0.5, -6.0, 0.4, 1.0, 0.2),
        (20, "T", "C", "MT-ATP6", "missense_variant", 0.3, -4.0, 0.2, 1.0, 0.2),          # one base, two genes
        (30, "G", "A", "MT-CO1", "synonymous_variant", np.nan, np.nan, np.nan, 2.0, 0.3),
        (40, "C", "T", "MT-CYB", "missense_variant", np.nan, -3.0, 0.3, 0.5, 0.1)]         # AlphaMissense gives none


class TestReadMitimpact:
    def test_missense_rows_with_the_13_names_and_the_signs(self, tmp_path):
        m = vc.read_mitimpact(mitimpact_csv(tmp_path, ROWS), digest=None)
        assert list(m.gene) == ["ND1", "ND1", "ATP8", "ATP6", "CYTB"] and set(m.gene) <= set(cybrid.ALL_GENES)
        first = m.iloc[0]
        assert (first.alphamissense, first.esm1b, first.apogee2, first.phylop_100v, first.mlc) == (0.9, 12.0, 0.8, 5.0, 0.7)
        assert np.isnan(m.alphamissense.iloc[4]) and m.esm1b.iloc[4] == 3.0

    def test_another_file_than_the_pinned_one_is_refused(self, tmp_path):
        with pytest.raises(SystemExit):
            vc.read_mitimpact(mitimpact_csv(tmp_path, ROWS))

    def test_a_change_twice_for_one_gene_and_an_unknown_gene_are_refused(self, tmp_path):
        with pytest.raises(ValueError, match="twice"):
            vc.read_mitimpact(mitimpact_csv(tmp_path, ROWS + [ROWS[0]]), digest=None)
        with pytest.raises(ValueError, match="not one of the 13"):
            vc.read_mitimpact(mitimpact_csv(tmp_path, [(5, "A", "C", "MT-XYZ", "missense_variant", 0.5, -1.0, 0.5, 1.0, 0.1)]), digest=None)


def variants():
    return pd.DataFrame({"pos": [10, 10, 20, 20, 40, 50], "ref": list("AATTCG"), "alt": list("CGCCTA"),
                         "gene": ["ND1", "ND1", "ATP8", "ATP6", "CYTB", "ND2"],
                         "label": ["harmful", "harmless", "harmful", "harmless", "harmless", "harmless"],
                         "n_aa": [1, 4, 2, 3, 5, 2], "stab_fold": [2.0, 0.1, 1.0, 0.5, 0.2, np.nan]})


KRON = pd.DataFrame({"pos": [10, 10, 20, 40, 60], "ref": list("AATCA"), "alt": list("CGCTG"), "kron": [5.0, 1.0, 3.0, -1.0, 9.0]})


class TestJoinedAndTable:
    def joined(self, tmp_path):
        return vc.joined(variants(), vc.read_mitimpact(mitimpact_csv(tmp_path, ROWS), digest=None), KRON)

    def test_the_join_is_by_position_bases_and_gene_and_a_missing_score_stays_empty(self, tmp_path):
        v = self.joined(tmp_path)
        assert len(v) == 6 and list(v.alphamissense[:4]) == [0.9, 0.1, 0.5, 0.3]          # ATP8 and ATP6 each get their own
        assert np.isnan(v.alphamissense[4]) and v.esm1b[4] == 3.0                         # CYTB: no AlphaMissense
        assert v.loc[5, ["alphamissense", "esm1b", "kron"]].isna().all()                  # ND2 at 50: in no table
        assert list(v.kron[:5]) == [5.0, 1.0, 3.0, 3.0, -1.0]                             # by base, so both genes at 20
        assert list(v[vc.BASELINE]) == [-1, -4, -2, -3, -5, -2] and list(v.harmful) == [True, False, True, False, False, False]

    def test_the_areas_are_on_the_variants_a_score_covers(self, tmp_path):
        t = vc.table(self.joined(tmp_path)).set_index("score")
        # AlphaMissense covers four: harmful 0.9 and 0.5 above harmless 0.1 and 0.3; the baseline there (-1, -2 above -4, -3) too
        assert (t.n_harmful["alphamissense"], t.n_harmless["alphamissense"], t.not_covered["alphamissense"]) == (2, 2, 2)
        assert t.area["alphamissense"] == 1.0 and t.baseline_area["alphamissense"] == 1.0 and t.difference["alphamissense"] == 0.0
        # ESM1b covers five; with its sign turned the harmful 12 and 6 are above 2, 4 and 3
        assert (t.n_harmful["esm1b"], t.n_harmless["esm1b"]) == (2, 3) and t.area["esm1b"] == 1.0
        assert t.low["esm1b"] <= t.difference["esm1b"] <= t.high["esm1b"]
        assert set(t.index) == set(vc.PREDICTORS) | set(vc.OURS)

    def test_a_score_upside_down_reads_below_one_half(self, tmp_path):
        v = self.joined(tmp_path)
        v["esm1b"] = -v.esm1b
        assert vc.table(v).set_index("score").area["esm1b"] == 0.0

    def test_one_kind_only_gets_no_area(self, tmp_path):
        v = self.joined(tmp_path)
        t = vc.table(v[~v.harmful]).set_index("score")
        assert t.area.isna().all() and (t.n_harmful == 0).all()
