"""The complex IV interface test, the structural-alignment check and Coevol's reader, on small real cases.

Expected values are worked out by hand in each test (or read from real program output kept in tests/data:
US-align on the human ND4 and ND5 chains, and the first Coevol chain of the published run)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines.longevity import coevol
from mitoage_pipelines.mitoport import burial, contacts, iface
from mitoage_pipelines.structure import fold, msa, within

DATA = Path(__file__).parent / "data"


class TestContactMap:
    def test_entity_descriptions_of_5z62_name_the_mitochondrial_subunits(self):
        # the descriptions as PDB 5Z62 gives them; a complex I accessory subunit must not be taken for COX1
        assert contacts.mt_gene("Cytochrome c oxidase subunit 1") == "COX1"
        assert contacts.mt_gene("Cytochrome c oxidase subunit 2") == "COX2"
        assert contacts.mt_gene("Cytochrome c oxidase subunit 3") == "COX3"
        assert contacts.mt_gene("Cytochrome c oxidase subunit 7C, mitochondrial") is None
        assert contacts.mt_gene("NADH dehydrogenase [ubiquinone] 1 alpha subcomplex subunit 1") is None
        assert burial.gene_of("Cytochrome c oxidase subunit 4 isoform 1, mitochondrial") == "COX4I1"
        assert burial.gene_of("Cytochrome c oxidase subunit NDUFA4") == "NDUFA4"
        assert burial.gene_of("Cytochrome c oxidase subunit 1") is None

    def test_column_of_position_skips_gaps(self):
        assert burial.column_of_position("-MA--K") == {1: 1, 2: 2, 3: 5}
        assert burial.column_of_position("") == {}

    def test_tree_tips_reads_slugs_with_hyphens(self):
        assert burial.tree_tips("((homo-sapiens:1.5,pan-troglodytes:1.5):2,mus-musculus:3.5);") == \
            {"homo-sapiens", "pan-troglodytes", "mus-musculus"}


class TestMatchedControls:
    def frame(self):
        # gene of six residues: two contacts (0.40, 0.10) and four candidates; the last has no column
        return pd.DataFrame({"n_contacts": [2, 0, 1, 0, 0, 0],
                             "rasa_alone": [0.40, 0.38, 0.10, 0.12, 0.90, 0.11],
                             "aln_col":    [5, 6, 7, 8, 9, 0]})

    def test_nearest_by_isolated_accessibility_without_replacement(self):
        con, taken = burial.match_controls(self.frame())
        # contacts in order of accessibility: 0.10 (index 2) first, then 0.40 (index 0)
        assert con == [2, 0]
        # 0.10 takes 0.12 (index 3; index 5 at 0.11 has no column); 0.40 takes 0.38 (index 1)
        assert taken == [3, 1]
        assert set(taken).isdisjoint(con)

    def test_an_exhausted_pool_leaves_contacts_unmatched(self):
        f = pd.DataFrame({"n_contacts": [1, 1, 0], "rasa_alone": [0.2, 0.3, 0.25], "aln_col": [1, 2, 3]})
        con, taken = burial.match_controls(f)
        assert con == [0, 1] and taken == [2]


class TestDistances:
    ALN = {"a": "MKVL-A", "b": "MKIL-A", "c": "MRIV-S", "d": "MK----"}

    def test_hamd_counts_differences_over_shared_columns_and_needs_five(self):
        D = iface.hamd(self.ALN, ["a", "b", "c", "d"], [0, 1, 2, 3, 4, 5])
        assert D[0, 1] == pytest.approx(1 / 5)          # a/b: V-I over the five shared columns
        assert D[0, 2] == pytest.approx(4 / 5)          # a/c: K-R, V-I, L-V, A-S
        assert np.isnan(D[0, 3])                        # a/d share two columns: too few
        assert np.allclose(np.diag(D)[:3], 0)

    def test_mismatch_is_against_the_modal_residue_and_ignores_gaps(self):
        share, filled = iface.mismatch(self.ALN, ["a", "b", "c", "d"], [1, 2, 5])
        # modes: column 1 K (3 of 4), column 2 I (2 of 3 filled; I beats V), column 5 A (2 of 3)
        assert share.tolist() == pytest.approx([1 / 3, 0, 2 / 3, 0])
        assert filled.tolist() == pytest.approx([1, 1, 1, 1 / 3])

    def test_cosine_distance_marks_unusable_rows(self):
        A = np.array([[1, 0], [0, 1], [np.nan, np.nan]], dtype=np.float32)
        D = iface.cosd(A)
        assert D[0, 1] == pytest.approx(1) and D[0, 0] == pytest.approx(0)
        assert np.isnan(D[2]).all() and np.isnan(D[:, 2]).all()

    def test_pcorr_removes_what_the_covariate_explains(self):
        rng = np.random.default_rng(1)
        T = rng.permutation(400) * 0.25        # evenly spaced, so its ranks are linear in it
        x = T + rng.normal(0, 1, 400); y = T + rng.normal(0, 1, 400)      # related only through T
        assert np.corrcoef(x, y)[0, 1] > 0.99
        assert abs(iface.pcorr(x, y, T)) < 0.15
        assert np.isnan(iface.pcorr(x[:150], y[:150], T[:150]))          # under 200 pairs


class TestPhylogeneticRegression:
    def test_gls_with_independent_tips_is_ordinary_least_squares(self):
        rng = np.random.default_rng(2)
        n = 50; X = np.column_stack([np.ones(n), rng.normal(size=n)]); y = 2 + 3 * X[:, 1] + rng.normal(0, 0.1, n)
        b, se, df = iface.gls(X, y, np.eye(n), lam=0.5)       # lambda scales off-diagonals only: none here
        ols, *_ = np.linalg.lstsq(X, y, rcond=None)
        assert b == pytest.approx(ols) and df == n - 2
        resid = y - X @ ols
        assert se[1] == pytest.approx(np.sqrt((resid @ resid / (n - 2)) * np.linalg.inv(X.T @ X)[1, 1]))

    def test_lambda_shrinks_shared_history_only(self):
        C = np.array([[2.0, 1.5], [1.5, 2.0]])
        assert iface.V(0.0, C).tolist() == [[2.0, 0.0], [0.0, 2.0]]
        assert iface.V(1.0, C).tolist() == C.tolist()


class TestStructuralCheck:
    def test_residue_pairs_follow_columns_not_positions(self):
        aln = {"a": "AB-C", "b": "A-DC"}
        assert msa.residue_pairs(aln, "a", "b") == {(0, 0), (2, 2)}

    def test_core_needs_fill_and_a_dominant_residue(self):
        rows = ["MKAW", "MKAW", "MRAW", "M-SW", "MKAW", "MKA-", "MKA-", "MKA-", "MKA-", "MKA-"]
        # column 0 all M; column 1 filled 9/10 with K 8/9 (below 0.9); column 2 A 9/10;
        # column 3 all W where filled, but filled in 5 of 10 rows only
        assert msa.core_columns(rows) == [0, 2]

    def test_foldmason_row_names_back_to_slugs(self):
        assert msa.foldmason_names({"acinonyx-jubatus.pdb_A": "M", "homo-sapiens.pdb": "K"}) == \
            {"acinonyx-jubatus": "M", "homo-sapiens": "K"}

    def test_usalign_report_of_nd4_on_nd5(self):
        r = within.parse_usalign((DATA / "usalign_ND4_ND5.txt").read_text())
        # 459 and 603 residues: normalised by the shorter (ND4) 0.89829, by the longer (ND5) 0.69351
        assert r == {"tm_short": 0.89829, "tm_long": 0.69351, "aligned": 447, "rmsd": 2.46}
        assert within.parse_usalign("no alignment here") is None

    def test_mean_plddt_reads_ca_b_factors(self):
        line = "ATOM      2  CA  MET A   1      11.104   6.134  -6.504  1.00 {:5.2f}           C"
        other = "ATOM      1  N   MET A   1      11.104   6.134  -6.504  1.00 10.00           N"
        text = "\n".join([line.format(90.0), other, line.format(80.0)])
        assert fold.mean_plddt(text) == pytest.approx(85.0)
        assert np.isnan(fold.mean_plddt(""))


class TestCoevolReader:
    def test_published_chain_one(self):
        b = coevol.parse_cov((DATA / "coevol_chain1.cov").read_text())
        # order dS, dN/dS, Longevity, Mass; the values the hypotheses document quotes as provisional
        assert b["correlation coefficients"][2][0] == pytest.approx(-0.518)        # dS-longevity
        assert b["correlation coefficients"][1][2] == pytest.approx(0.449)         # dN/dS-longevity
        assert b["posterior probs"][1][2] == pytest.approx(0.99)
        assert b["partial correlation coefficients"][1][2] == pytest.approx(0.283)
        assert b["partial posterior probs"][0][2] == pytest.approx(0.067)
        assert b["posterior probs"][0][0] is None
        for k in ("covariances", "correlation coefficients", "posterior probs", "precisions",
                  "partial correlation coefficients", "partial posterior probs"):
            assert len(b[k]) == 4 and all(len(r) == 4 for r in b[k])
