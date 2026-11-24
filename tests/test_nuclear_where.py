"""`mitoport nuclear-where` on frames small enough to count by hand."""
import numpy as np
import pandas as pd

from mitoage_pipelines.mitoport import nuclear_where as nw

CHAINS = pd.DataFrame({"chain": ["A", "B"], "gene": ["NUCA", "NUCB"]})
CONTACTS = pd.DataFrame({"mt_gene": ["G", "G", "G", "G"], "human_pos": [1, 1, 2, 3], "partner_chain": ["A", "B", "A", "A"],
                         "partner_is_mt": [0, 0, 0, 1], "partner_resnum": [10, 20, 11, 99], "partner_res": list("LKVA")})
KEY = ["host", "donor", "complex"]


def diff(rows):
    return pd.DataFrame([{"host": "h", "donor": d, "complex": "X", "gene": "G", "human_pos": p} for d, p in rows])


class TestFacing:
    def test_a_difference_faces_the_nuclear_residues_that_touch_it_and_not_a_mitochondrial_one(self):
        f = nw.facing(diff([("d1", 1), ("d1", 3), ("d2", 2)]), CONTACTS, CHAINS)
        got = sorted(zip(f.donor, f.nuclear_gene, f.partner_resnum))
        assert got == [("d1", "NUCA", 10), ("d1", "NUCB", 20), ("d2", "NUCA", 11)]        # position 3 touches only a mitochondrial chain

    def test_a_residue_is_counted_once_per_pair(self):
        d = pd.concat([diff([("d1", 1)]), diff([("d1", 1)])])
        assert len(nw.facing(d, CONTACTS, CHAINS)) == 2


class TestSubunitsAndResidues:
    lab = pd.DataFrame({"host": "h", "donor": ["d1", "d2", "d3"], "complex": "X", "status": ["reduced", "not_reduced", "not_reduced"],
                        "pct": [50.0, 95.0, 100.0]})
    d = diff([("d1", 1), ("d1", 2), ("d2", 2), ("d3", 3)])

    def test_the_count_separates_reduced_from_not(self):
        face = nw.facing(self.d, CONTACTS, CHAINS)
        s = nw.subunits(self.lab, self.d, face).set_index("nuclear_gene")
        assert s.loc["NUCA", "mean_facing_reduced"] == 2 and s.loc["NUCA", "mean_facing_not_reduced"] == 0.5
        assert s.loc["NUCA", "auc"] == 1.0 and s.loc["NUCB", "mean_facing_reduced"] == 1
        assert s.loc["NUCA", "pairs_reduced"] == 1 and s.loc["NUCA", "pairs_not_reduced"] == 2

    def test_residue_shares(self):
        r = nw.residues(self.lab, nw.facing(self.d, CONTACTS, CHAINS)).set_index(["nuclear_gene", "partner_resnum"])
        assert r.loc[("NUCA", 10), "share_reduced"] == 1 and r.loc[("NUCA", 10), "share_not_reduced"] == 0
        assert r.loc[("NUCA", 11), "share_reduced"] == 1 and r.loc[("NUCA", 11), "share_not_reduced"] == 0.5
        assert nw.residues(self.lab, nw.facing(self.d, CONTACTS, CHAINS)).iloc[0]["difference"] == 1.0


class TestLabels:
    def test_disagreeing_papers_are_left_out_and_100_percent_is_not_reduced(self):
        m = pd.DataFrame({"host": "h", "donor": ["a", "b", "b", "c"], "complex": "X", "significant": ["yes", "yes", "no", "yes"],
                          "pct": [50.0, 60.0, 90.0, 100.0]})
        lab = nw.labelled(m).set_index("donor")
        assert list(lab.index) == ["a", "c"] and lab.loc["a", "status"] == "reduced" and lab.loc["c", "status"] == "not_reduced"


class TestOrder:
    def test_residues_are_ranked_by_the_difference_of_shares_not_by_how_often_a_reduced_pair_faces_them(self):
        lab = pd.DataFrame({"host": "h", "donor": ["d1", "d2", "d3", "d4"], "complex": "X",
                            "status": ["reduced", "reduced", "reduced", "not_reduced"], "pct": [50.0, 55.0, 60.0, 98.0]})
        d = diff([("d1", 1), ("d2", 1), ("d4", 1), ("d3", 2)])
        r = nw.residues(lab, nw.facing(d, CONTACTS, CHAINS))
        # residue 10 faces 2 of 3 reduced and 1 of 1 not reduced (-1/3); residue 11 faces 1 of 3 reduced and 0 of 1 (+1/3)
        assert list(r.partner_resnum[r.nuclear_gene == "NUCA"]) == [11, 10]
        assert np.allclose(r.difference[r.nuclear_gene == "NUCA"], [1 / 3, -1 / 3])
