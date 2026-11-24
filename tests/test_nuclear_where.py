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


class TestResidues:
    lab = pd.DataFrame({"host": "h", "donor": ["d1", "d2", "d3"], "complex": "X", "status": ["reduced", "not_reduced", "not_reduced"],
                        "pct": [50.0, 95.0, 100.0]})
    d = diff([("d1", 1), ("d1", 2), ("d2", 2), ("d3", 3)])

    def test_residue_shares(self):
        r = nw.residues(self.lab, nw.facing(self.d, CONTACTS, CHAINS)).set_index(["nuclear_gene", "partner_resnum"])
        assert r.loc[("NUCA", 10), "share_reduced"] == 1 and r.loc[("NUCA", 10), "share_not_reduced"] == 0
        assert r.loc[("NUCA", 11), "share_reduced"] == 1 and r.loc[("NUCA", 11), "share_not_reduced"] == 0.5
        assert nw.residues(self.lab, nw.facing(self.d, CONTACTS, CHAINS)).iloc[0]["difference"] == 1.0


class TestMasksAndCounts:
    def test_masks_count_distinct_residues_of_each_subunit(self):
        m = nw.site_masks(CONTACTS, CHAINS)
        assert set(m) == {("G", 1), ("G", 2)}                               # position 3 touches only a mitochondrial chain
        assert nw.pair_counts([("G", 1)], m) == {"NUCA": 1, "NUCB": 1}
        assert nw.pair_counts([("G", 1), ("G", 2)], m) == {"NUCA": 2, "NUCB": 1}
        assert nw.pair_counts([("G", 1), ("G", 1)], m) == {"NUCA": 1, "NUCB": 1}   # a residue is counted once
        assert nw.pair_counts([("G", 9)], m) == {}

    def test_a_residue_listed_twice_for_one_site_is_one_bit(self):
        c = pd.concat([CONTACTS, CONTACTS.iloc[[0]]])
        m = nw.site_masks(c, CHAINS)
        assert nw.pair_counts([("G", 1)], m) == {"NUCA": 1, "NUCB": 1}
        # the repeated row must not spill into the bit of the subunit's next residue (the one site 2 touches)
        assert nw.pair_counts([("G", 1), ("G", 2)], m) == {"NUCA": 2, "NUCB": 1}

    def test_two_sites_touching_the_same_residue_count_it_once(self):
        c = pd.DataFrame({"mt_gene": ["G", "G"], "human_pos": [1, 2], "partner_chain": ["A", "A"], "partner_is_mt": [0, 0],
                          "partner_resnum": [10, 10], "partner_res": ["L", "L"]})
        assert nw.pair_counts([("G", 1), ("G", 2)], nw.site_masks(c, CHAINS)) == {"NUCA": 1}

    def test_area_columns_equal_the_area_of_cybrid(self):
        from mitoage_pipelines.mitoport import cybrid
        x = np.array([[3, 0], [1, 2], [2, 2], [0, 2], [1, 0]], float)
        red = np.array([True, True, False, False, False])
        got = nw.auc_columns(x, red)
        for k in range(2):
            assert got[k] == cybrid.auc(pd.Series(x[:, k]), pd.Series(red))
        assert np.isnan(nw.auc_columns(x, np.zeros(5, bool))).all()


class TestNull:
    def test_draws_keep_the_number_of_sites_per_protein_and_take_them_from_the_universe(self):
        sites = {("h", "d", "X"): [("G", 1), ("G", 2), ("H", 5)]}
        universe = {"G": np.array([1, 2, 3, 4, 5, 6]), "H": np.array([5, 7, 9])}
        rng = np.random.default_rng(0)
        for _ in range(20):
            d = nw.draw_sites(sites, universe, rng)[("h", "d", "X")]
            assert sorted(g for g, _ in d) == ["G", "G", "H"] and len({p for g, p in d if g == "G"}) == 2
            assert all(p in universe[g] for g, p in d)
        assert nw.draw_sites(sites, universe, np.random.default_rng(1)) == nw.draw_sites(sites, universe, np.random.default_rng(1))


class TestSubunits:
    contacts = pd.DataFrame({"mt_gene": ["G"] * 4, "human_pos": [1, 2, 3, 4], "partner_chain": ["A", "A", "B", "B"],
                             "partner_is_mt": [0] * 4, "partner_resnum": [10, 11, 20, 21], "partner_res": list("LKVA")})
    chains = pd.DataFrame({"chain": ["A", "B"], "gene": ["NUCA", "NUCB"]})
    lab = pd.DataFrame({"host": "h", "donor": [f"d{i}" for i in range(1, 8)], "complex": "X", "pct": [40, 45, 50, 90, 95, 100.0, 100.0],
                        "status": ["reduced"] * 3 + ["not_reduced"] * 4})              # d7 differs nowhere
    # the reduced donors differ at sites 1 and 2 (the surface of NUCA); the others at 3 (NUCB)
    d = diff([("d1", 1), ("d1", 2), ("d2", 1), ("d2", 2), ("d3", 1), ("d3", 2), ("d4", 3), ("d5", 3), ("d6", 3)])

    def run(self, draws=200, lab=None):
        masks = {"X": nw.site_masks(self.contacts, self.chains)}
        universe = {"X": {"G": np.array([1, 2, 3, 4])}}
        return nw.subunits(self.lab if lab is None else lab, self.d, masks, universe, draws, seed=0)

    def test_the_subunit_facing_the_reduced_donors_changes_has_the_top_area(self):
        s, expo = self.run()
        s = s.set_index("nuclear_gene")
        assert s.loc["NUCA", "auc"] == 1.0 and s.loc["NUCB", "auc"] == 0.125      # d7 differs nowhere: three ties of twelve
        assert s.loc["NUCA", "mean_facing_reduced"] == 2 and s.loc["NUCA", "mean_facing_not_reduced"] == 0
        assert list(expo.loc[expo.donor == "d1", ["NUCA", "NUCB"]].iloc[0]) == [2, 0]
        assert list(expo.loc[expo.donor == "d4", ["NUCA", "NUCB"]].iloc[0]) == [0, 1]

    def test_the_null_is_the_interfaces_own_expectation(self):
        s, _ = self.run()
        s = s.set_index("nuclear_gene")
        assert 0.0 < s.loc["NUCA", "auc_null_mean"] < 1.0 and 0.0 <= s.loc["NUCA", "p_subunit"] <= 1.0
        assert s.loc["NUCA", "p_family"] == s.loc["NUCB", "p_family"] and s.loc["NUCA", "p_family"] > 0
        # only 4 positions and 3 sites a pair: the best area is reached by chance in many draws, so nothing is named
        assert not s.named.any()

    def test_leaving_a_donor_out_gives_a_range(self):
        # d3 is reduced but differs only at the surface of NUCB: NUCA's counts are [2, 2, 0] for the reduced donors
        d = diff([("d1", 1), ("d1", 2), ("d2", 1), ("d2", 2), ("d3", 3), ("d4", 3), ("d5", 3), ("d6", 3)])
        masks = {"X": nw.site_masks(self.contacts, self.chains)}
        s, _ = nw.subunits(self.lab, d, masks, {"X": {"G": np.array([1, 2, 3, 4])}}, 20, seed=0)
        a = s.set_index("nuclear_gene").loc["NUCA"]
        assert np.isclose(a.auc, (2 * 4 + 4 * 0.5) / 12)          # two reduced donors win all four comparisons, d3 ties them
        assert a.auc_leave_donor_min < a.auc < a.auc_leave_donor_max          # without d1 or d2 lower, without d3 higher
        assert a.auc_leave_donor_max == 1.0

    def test_p_values_count_the_observed_draw_too(self):
        s, _ = self.run(draws=200)
        for col in ("p_subunit", "p_family"):
            v = s[col].to_numpy() * 201
            assert np.allclose(v, np.round(v)) and (v >= 1 - 1e-9).all()

    def test_identical_columns_share_a_group(self):
        c = pd.concat([self.contacts, pd.DataFrame({"mt_gene": ["G"], "human_pos": [1], "partner_chain": ["C"], "partner_is_mt": [0],
                                                    "partner_resnum": [30], "partner_res": ["M"]}),
                       pd.DataFrame({"mt_gene": ["G"], "human_pos": [2], "partner_chain": ["C"], "partner_is_mt": [0],
                                     "partner_resnum": [31], "partner_res": ["M"]})])
        chains = pd.DataFrame({"chain": ["A", "B", "C"], "gene": ["NUCA", "NUCB", "NUCC"]})
        masks = {"X": nw.site_masks(c, chains)}
        s, _ = nw.subunits(self.lab, self.d, masks, {"X": {"G": np.array([1, 2, 3, 4])}}, 50, seed=0)
        s = s.set_index("nuclear_gene")
        assert s.loc["NUCA", "profile_group"] == s.loc["NUCC", "profile_group"] != s.loc["NUCB", "profile_group"]
        assert s.loc["NUCA", "profile_group_size"] == 2 and s.loc["NUCB", "profile_group_size"] == 1

    def test_a_subunit_the_reduced_donors_change_more_than_chance_is_named_and_its_p_is_family_wise(self):
        # 40 positions: 1-10 touch NUCA, 11-20 touch NUCB, 21-40 touch nothing; the reduced donors differ at five sites of
        # NUCA's surface, the others at five of NUCB's
        c = pd.DataFrame({"mt_gene": "G", "human_pos": list(range(1, 21)), "partner_chain": ["A"] * 10 + ["B"] * 10,
                          "partner_is_mt": 0, "partner_resnum": list(range(1, 21)), "partner_res": "L"})
        d = diff([(f"d{i}", p) for i in (1, 2, 3) for p in (1, 2, 3, 4, 5)] + [(f"d{i}", p) for i in (4, 5, 6, 7) for p in (11, 12, 13, 14, 15)])
        masks = {"X": nw.site_masks(c, self.chains)}
        s, _ = nw.subunits(self.lab, d, masks, {"X": {"G": np.arange(1, 41)}}, 500, seed=0)
        s = s.set_index("nuclear_gene")
        assert s.loc["NUCA", "auc"] == 1.0 and s.loc["NUCB", "auc"] == 0.0
        assert s.loc["NUCA", "p_family"] < 0.05 and s.loc["NUCA", "named"] and not s.loc["NUCB", "named"]
        assert s.loc["NUCA", "p_family"] >= s.loc["NUCA", "p_subunit"]

    def test_the_same_seed_gives_the_same_table(self):
        a, _ = self.run(); b, _ = self.run()
        assert a.equals(b)


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

    def test_disagreeing_papers_can_be_coded_reduced(self):
        m = pd.DataFrame({"host": "h", "donor": ["b", "b"], "complex": "X", "significant": ["yes", "no"], "pct": [60.0, 90.0]})
        assert len(nw.labelled(m)) == 0 and nw.labelled(m, "reduced").status.iloc[0] == "reduced"

    def test_rows_of_other_complexes_are_dropped_when_asked(self):
        m = pd.DataFrame({"host": "h", "donor": ["a", "a"], "complex": ["X", "all"], "significant": ["yes", "yes"], "pct": [50.0, 50.0]})
        assert list(nw.labelled(m, complexes=("X",))["complex"]) == ["X"] and list(nw.labelled(m)["complex"]) == ["X", "all"]
