"""H16, the simulated cybrids: gene names, placing a protein on human positions, the counts, the agreement of
sources, the ROC area, the shipped measurement table, and the complex IV contact map against `mitoport contacts`.

Expected values are worked out by hand in each test; the structure test compares two code paths on the real
file and runs only where `download pdb 5Z62` has put it."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import complexes, contacts, cybrid


class TestNames:
    def test_spellings_of_the_13_genes(self):
        assert cybrid.gene_name("COX1", None) == "COX1"
        assert cybrid.gene_name("COI", None) == "COX1"
        assert cybrid.gene_name("mt-Cytb", None) == "CYTB"
        assert cybrid.gene_name("cytb", None) == "CYTB"
        assert cybrid.gene_name("NAD4L", None) == "ND4L"
        assert cybrid.gene_name(None, "cytochrome b") == "CYTB"      # AY193770 has a product and no gene
        assert cybrid.gene_name("mt-Rnr1", None) is None
        assert cybrid.gene_name(None, None) is None


class TestOnHuman:
    def test_a_substitution_and_an_internal_deletion(self):
        human = "MKVLAWYRTE"
        other = "MKVLAYRSE"                       # W deleted, T -> S
        got = cybrid.on_human(human, other, partial=False)
        assert "".join(got) == "MKVLA-YRSE"

    def test_the_ends_of_a_partial_sequence_are_unknown_not_deleted(self):
        human = "MKVLAWYRTEQ"
        other = "VLAWYRTE"                        # starts at human position 3, stops before the last two
        assert "".join(cybrid.on_human(human, other, partial=True)) == "???LAWYRT??"
        assert "".join(cybrid.on_human(human, other, partial=False)) == "--VLAWYRTE-"

    def test_the_first_residue_of_a_partial_sequence_is_not_a_substitution(self):
        human = "MPMANLLLLI"
        other = "TPMANLLLLI"                      # a record partial at its start: T is not the protein's first residue
        assert cybrid.on_human(human, other, partial=True)[0] == "?"
        assert cybrid.on_human(human, other, partial=False)[0] == "T"


class TestDonorThroughHost:
    # the host and the donor carry the same three extra residues that human lacks; aligned to human one by one the
    # two gaps can fall in different places, aligned to each other they are the same residues
    human = "MKVLAWYRTEQNDPSG"
    host = "MKVLAWAAAYRTEQNDPSG"
    donor = "MKVLAWAAAYRSEQNDPSG"                  # one real difference: T -> S

    def test_the_donor_is_read_against_the_host(self):
        hs = cybrid.on_human(self.human, self.host, False)
        ds = cybrid.donor_through_host(self.human, self.host, self.donor, False)
        assert sum(a != b for a, b in zip(hs, ds)) == 1
        assert [(a, b) for a, b in zip(hs, ds) if a != b] == [("T", "S")]

    def test_a_gap_in_the_donor_alone_is_a_difference_once(self):
        hs = cybrid.on_human("MKVLAWYRTE", "MKVLAWYRTE", False)
        ds = cybrid.donor_through_host("MKVLAWYRTE", "MKVLAWYRTE", "MKVLAYRTE", False)
        assert [(a, b) for a, b in zip(hs, ds) if a != b] == [("W", "-")]

    def test_a_shifted_gap_is_not_counted(self):
        # host and donor are equal and both lack one of human's two tryptophans: whichever one each alignment to
        # human drops, the pair has no difference
        human, same = "MKVLAWWYRTE", "MKVLAWYRTE"
        hs = cybrid.on_human(human, same, False); ds = cybrid.donor_through_host(human, same, same, False)
        assert hs == ds

    # ATP8 of the three reference records (NC_012920.1, NC_005089.1, NC_025268.1): the case the review found
    ATP8 = {"human": "MPQLNTTVWPTMITPMLLTLFLITQLKMLNTNYHLPPSPKPMKMKNYNKPWEPKWTKICSLHSLPPQS",
            "mouse": "MPQLDTSTWFITIISSMITLFILFQLKVSSQTFPLAPSPKSLTTMKVKTPWELKWTKIYLPHSLPQQ",
            "caroli": "MPQLDTSTWFITIISSMITLFILFQLKISSQIFPLPPSPKSLTTMKVKNPWESKWTKIYLPHSLPQQ"}

    def test_mouse_and_caroli_atp8_differ_at_five_residues(self):
        h, m, c = self.ATP8["human"], self.ATP8["mouse"], self.ATP8["caroli"]
        direct = sum(x != y for x, y in zip(m, c))
        assert len(m) == len(c) and direct == 5
        hs = cybrid.on_human(h, m, False)
        assert sum(x != y for x, y in zip(hs, cybrid.donor_through_host(h, m, c, False))) == direct
        # each aligned to human by itself, the same pair counts more: what the count did before
        assert sum(x != y for x, y in zip(hs, cybrid.on_human(h, c, False))) > direct

    def test_the_count_reads_the_donor_through_the_host(self):
        h, m, c = self.ATP8["human"], self.ATP8["mouse"], self.ATP8["caroli"]
        positions = pd.DataFrame({"gene": "ATP8", "human_pos": range(1, len(h) + 1), "resolved": 1, "contact_nuclear": 0})
        s = cybrid.pair_scores(positions, {"ATP8": h}, {"ATP8": (m, False, "m")}, {"ATP8": (c, False, "c")})
        assert s["d_all"] == sum(x != y for x, y in zip(m, c)) and s["n_unknown"] == 0

    def test_index_on_human(self):
        assert cybrid.index_on_human("MKVLAWYRTE", "MKVLAYRTE") == [0, 1, 2, 3, 4, None, 5, 6, 7, 8]


class TestPairScores:
    def positions(self):
        # gene X, five human positions: 1 and 2 touch a nuclear subunit, 3 and 4 do not, 5 is not resolved
        return pd.DataFrame({"gene": ["X"] * 5, "human_pos": [1, 2, 3, 4, 5],
                             "resolved": [1, 1, 1, 1, 0], "contact_nuclear": [1, 1, 0, 0, 0]})

    def test_counts_by_hand(self):
        human = {"X": "MKWLA"}
        host = {"X": ("MKWLA", False, "h")}
        donor = {"X": ("MRWLG", False, "d")}       # K->R at a contact (BLOSUM62 +2), A->G unresolved
        s = cybrid.pair_scores(self.positions(), human, host, donor)
        assert (s["n_contact"], s["d_contact"], s["d_contact_radical"]) == (2, 1, 0)
        assert (s["n_noncontact"], s["d_noncontact"]) == (2, 0)
        assert (s["n_all"], s["d_all"], s["n_unknown"]) == (5, 2, 0)
        assert s["f_contact"] == 0.5 and s["f_all"] == 0.4

    def test_a_radical_change_and_an_unknown_position(self):
        human = {"X": "MKWLA"}
        host = {"X": ("MKWLA", False, "h")}
        donor = {"X": ("MDWL", True, "d")}         # K->D at a contact (BLOSUM62 -1); partial, so position 5 and the
        #                                            record's own first and last residue (1 and 4) are unknown
        s = cybrid.pair_scores(self.positions(), human, host, donor)
        assert (s["d_contact"], s["d_contact_radical"]) == (1, 1)
        assert (s["n_all"], s["n_unknown"]) == (2, 3)

    def test_identical_proteins_score_zero(self):
        human = {"X": "MKWLA"}
        same = {"X": ("MKWLA", False, "x")}
        s = cybrid.pair_scores(self.positions(), human, same, same)
        assert s["d_contact"] == s["d_noncontact"] == s["d_all"] == 0


class TestStatusAndAuc:
    def rows(self, *pairs):
        return pd.DataFrame([{"significant": sig, "pct": pct} for sig, pct in pairs])

    def test_sources_that_agree_and_disagree(self):
        assert cybrid.status(self.rows(("yes", 46.0), ("yes", 51.5))) == "reduced"
        assert cybrid.status(self.rows(("no", np.nan))) == "not_reduced"
        assert cybrid.status(self.rows(("yes", 37.0), ("no", np.nan))) == "disagree"
        assert cybrid.status(self.rows(("yes", 120.0))) == "not_reduced"     # a significant increase is no defect

    def test_auc_by_hand(self):
        reduced = pd.Series([True, True, False, False])
        assert cybrid.auc(pd.Series([3, 4, 1, 2]), reduced) == 1.0
        assert cybrid.auc(pd.Series([1, 2, 3, 4]), reduced) == 0.0
        assert cybrid.auc(pd.Series([2, 2, 2, 1]), reduced) == 0.75          # 3 wins and 1 tie of 4 pairs
        assert np.isnan(cybrid.auc(pd.Series([1.0]), pd.Series([True])))


class TestMeasurementTable:
    def test_the_shipped_table(self):
        m = cybrid.read_measurements()
        assert set(m.host) | set(m.donor) <= set(cybrid.GENOMES)
        assert set(m.complex) <= set(cybrid.STRUCTURES) | {"all"}
        assert (m.source_doi.str.startswith("10.")).all() and (m.location != "").all()
        stated = m[m.pct.notna()]
        assert ((stated.pct_low <= stated.pct_high) & (stated.pct_low >= 0)).all()
        # a value is stated only for a measurement; a failed cybrid has none
        assert m[m.outcome == "no_respiring_cybrid"].pct.isna().all()
        assert set(m.significant) <= {"yes", "no", "not_stated"}


@pytest.mark.parametrize("complex_", ["I", "III", "IV", "V"])
def test_contacts_agree_with_mitoport_contacts(complex_, tmp_path):
    # two code paths over the same file: human positions from NCBI's human genome here, from the MitoAge
    # alignment in `mitoport contacts`; the contact residues must be the same set
    cx = complexes.get(complex_)
    if not cx.cif.is_file():
        pytest.skip(f"needs `mitoage download pdb {cx.pdb}`")
    if not (paths.CYBRID_GB / "records.gb").is_file():
        pytest.skip("needs `mitoage download cybrid`")
    if not (paths.BUILD / "alignments" / "COX1.aa.fasta.gz").is_file():
        pytest.skip("needs `mitoage download build` for mitoport contacts")
    records = cybrid.read_proteins(paths.CYBRID_GB / "records.gb")
    human = {g: p for g, (p, _) in records[cybrid.GENOMES["homo-sapiens"]].items()}
    mine = cybrid.contact_map(cx.cif, complex_, human)
    mine = {(g, p) for g, p, c in zip(mine.gene, mine.human_pos, mine.contact_nuclear) if c}
    table = pd.read_csv(contacts.main(cx.cif, paths.BUILD, tmp_path, cx), sep="\t")
    nuc = table[(table.partner_is_mt == 0) & table.human_pos.notna()]
    theirs = {(g, int(p)) for g, p in zip(nuc.mt_gene, nuc.human_pos)}
    assert mine == theirs and len(mine) > 0
