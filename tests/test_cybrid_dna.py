"""The nucleotide baselines beside H16 (exploratory): counting differences, carrying a place of the human rRNA gene
to another sequence, and, where `download cybrid` has put the records, the reading frames of humanin and MOTS-c on
the human record and a pair's rows.

Expected values are worked out by hand in each test, or read from the record itself."""
import pytest
from Bio.Seq import Seq

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid, cybrid_dna

GB = paths.CYBRID_GB / "records.gb"
needs_records = pytest.mark.skipif(not GB.is_file(), reason="needs `mitoage download cybrid`")


class TestDifferences:
    al = cybrid_dna.aligner()

    def test_the_same_sequence_has_none(self):
        assert cybrid_dna.differences("ACGTACGTAC", "ACGTACGTAC", self.al) == (0, 10)

    def test_a_mismatch_and_a_gap_each_count_once(self):
        a = "ACGTACGTACGTACGTACGT"
        assert cybrid_dna.differences(a, a[:9] + "T" + a[10:], self.al) == (1, 20)
        assert cybrid_dna.differences(a, a[:9] + a[10:], self.al) == (1, 20)

    def test_an_empty_side_differs_everywhere(self):
        assert cybrid_dna.differences("ACG", "", self.al) == (3, 3)


class TestCarried:
    al = cybrid_dna.aligner()
    human = "AAAAACCCCCGGGGGTTTTTAAAAACCCCC"

    def test_the_same_places_without_a_gap(self):
        other = self.human[:12] + "A" + self.human[13:]
        assert cybrid_dna.carried(self.human, other, 10, 15, self.al) == "GGAGG"

    def test_a_deletion_inside_shortens_it_and_one_outside_moves_it(self):
        inside = self.human[:12] + self.human[13:]
        assert cybrid_dna.carried(self.human, inside, 10, 15, self.al) == "GGGG"
        before = self.human[:2] + self.human[3:]
        assert cybrid_dna.carried(self.human, before, 10, 15, self.al) == "GGGGG"

    def test_an_insertion_inside_is_kept(self):
        other = self.human[:12] + "T" + self.human[12:]
        assert cybrid_dna.carried(self.human, other, 10, 15, self.al) == "GGTGGG"

    def test_an_insertion_right_after_the_last_place_is_not(self):
        other = self.human[:15] + "A" + self.human[15:]
        assert cybrid_dna.carried(self.human, other, 10, 15, self.al) == "GGGGG"
        assert cybrid_dna.carried(self.human, other, 15, 20, self.al) == "TTTTT"


class TestCdsOf:
    def test_the_cell_lines_own_record_where_one_is_listed_else_the_genome(self):
        (slug, gene), own = next(iter(cybrid.OWN_LINE.items()))
        features = {own: {"cds": {gene: "AAA"}}, cybrid.GENOMES[slug]: {"cds": {gene: "CCC", "COX1": "GGG"}}}
        assert cybrid_dna.cds_of(features, slug, gene) == "AAA"
        assert cybrid_dna.cds_of(features, slug, "COX1") == "GGG"


@pytest.fixture(scope="module")
def features():
    return cybrid_dna.read_features(GB)


@needs_records
class TestOnTheRecords:
    def test_the_frames_on_the_human_record_are_the_two_peptides(self, features):
        human = features[cybrid.GENOMES[cybrid.HUMAN]]
        f = cybrid_dna.frames(features, cybrid.HUMAN, human, cybrid_dna.aligner())
        assert str(Seq(f["humanin"]).translate(table=1)) == "MAPRGFSCLLLLTSEIDLPVKRRA*"
        assert str(Seq(f["motsc"]).translate(table=1)) == "MRWQEMGYIFYPRKLR*"

    def test_every_genome_has_two_rrna_genes_and_22_trna_genes(self, features):
        for acc in cybrid.GENOMES.values():
            assert set(features[acc]["rrna"]) == {"12S", "16S"}
            assert len(features[acc]["trna"]) == 22

    def test_a_species_against_itself_scores_nothing(self, features):
        for r in cybrid_dna.pair_rows(features, "mus-musculus", "mus-musculus"):
            assert all(r[s] == 0 for s in cybrid_dna.SCORES)

    def test_a_pair_is_the_same_both_ways_and_farther_donors_differ_more(self, features):
        ab = cybrid_dna.pair_rows(features, "homo-sapiens", "pan-troglodytes")
        ba = cybrid_dna.pair_rows(features, "pan-troglodytes", "homo-sapiens")
        far = cybrid_dna.pair_rows(features, "homo-sapiens", "pongo-abelii")
        assert [r["complex"] for r in ab] == list(cybrid.COMPLEX_GENES)
        for x, y, z in zip(ab, ba, far):
            assert x["nt_cds_d"] == y["nt_cds_d"] and x["nt_rrna_d"] == y["nt_rrna_d"]
            assert 0 < x["nt_cds_d"] < z["nt_cds_d"]
            assert 0 < x["nt_rrna_f"] < z["nt_rrna_f"] < 1
