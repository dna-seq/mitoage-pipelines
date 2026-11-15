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
        assert "".join(cybrid.on_human(human, other, partial=True)) == "??VLAWYRTE?"
        assert "".join(cybrid.on_human(human, other, partial=False)) == "--VLAWYRTE-"


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
        donor = {"X": ("MDWL", True, "d")}         # K->D at a contact (BLOSUM62 -1); partial, so position 5 unknown
        s = cybrid.pair_scores(self.positions(), human, host, donor)
        assert (s["d_contact"], s["d_contact_radical"]) == (1, 1)
        assert (s["n_all"], s["n_unknown"]) == (4, 1)

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
