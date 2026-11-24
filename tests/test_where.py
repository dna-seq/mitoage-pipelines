"""The exploratory tables of `mitoport where`, on frames small enough to count by hand."""
import numpy as np
import pandas as pd

from mitoage_pipelines.mitoport import where

SITES = pd.DataFrame({
    "complex": ["X"] * 6, "gene": ["G"] * 6, "human_pos": [1, 2, 3, 4, 5, 6], "human_aa": list("MKVLAW"),
    "n_aa": [1, 2, 3, 6, np.nan, 1], "rasa_complex": [0.05, 0.5, 0.1, np.nan, 0.3, 0.19],
    "conserved": [1, 1, 0, 0, np.nan, 1], "buried": [1, 0, 1, np.nan, 0, 1],
    "contact_nuclear": [1, 0, 0, 0, 0, 1], "contact_mt": [0, 0, 1, 0, 0, 0], "nuclear_chains": ["a", "", "", "", "", "b"]})


class TestVariantsWhere:
    v = pd.DataFrame({"id": list("abcde"), "gene": ["G"] * 5, "protein_pos": [1, 6, 4, 3, 5],
                      "label": ["harmful", "harmful", "harmful", "harmless", "harmless"], "n_aa": [1, 1, 6, 3, np.nan]})

    def test_shares_over_the_variants_whose_flag_is_known(self):
        t = where.variants_where(self.v, SITES).set_index(["what", "label"])
        assert tuple(t.loc[("conserved", "harmful"), ["known", "with"]]) == (3, 2)
        assert tuple(t.loc[("conserved", "harmless"), ["known", "with"]]) == (1, 0)          # position 5 is unknown
        assert tuple(t.loc[("buried", "harmful"), ["known", "with"]]) == (2, 2)              # position 4 has no residue
        assert tuple(t.loc[("contact_nuclear", "harmful"), ["known", "with"]]) == (3, 2)
        assert tuple(t.loc[("amino acids at the site: 1", "harmful"), ["known", "with"]]) == (3, 2)
        assert tuple(t.loc[("amino acids at the site: 5 or more", "harmful"), ["known", "with"]]) == (3, 1)
        assert tuple(t.loc[("amino acids at the site: 3 to 4", "harmless"), ["known", "with"]]) == (1, 1)


class TestCybrid:
    human = {"G": "MKVLAW"}
    prot = {"h": {"G": ("MKVLAW", False, "x")}, "d": {"G": ("MRILGF", False, "y")}}       # K2R, V3I, A5G, W6F

    def test_one_row_per_difference_with_its_site(self):
        d = where.cybrid_sites(SITES, self.human, self.prot, [("h", "d")])
        assert list(zip(d.human_pos, d.host_aa, d.donor_aa)) == [(2, "K", "R"), (3, "V", "I"), (5, "A", "G"), (6, "W", "F")]
        assert list(d.blosum) == [2.0, 3.0, 0.0, 1.0] and list(d.conserved[[0, 1, 3]]) == [1, 0, 1] and np.isnan(d.conserved[2])
        same = where.cybrid_sites(SITES, self.human, {"h": self.prot["h"], "d": self.prot["h"]}, [("h", "d")])
        assert len(same) == 0

    def test_the_summary_counts_conserved_sites(self):
        d = where.cybrid_sites(SITES, self.human, self.prot, [("h", "d")])
        meas = pd.DataFrame({"host": ["h"], "donor": ["d"], "complex": ["X"], "outcome": ["measured"], "pct": [60.0], "significant": ["yes"]})
        s = where.cybrid_summary(d, meas).iloc[0]
        assert (s.differences, s.at_conserved, s.at_conserved_buried, s.at_conserved_contact_nuclear, s.at_contact_nuclear) == (4, 2, 1, 1, 1)          # position 3 is buried and not conserved
        assert s.at_conserved_radical == 0 and s.pct == 60.0
