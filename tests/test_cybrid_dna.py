"""The nucleotide baselines beside H16 (exploratory): counting differences, carrying a place of the human rRNA gene
to another sequence, and, where `download cybrid` has put the records, the reading frames of humanin and MOTS-c on
the human record and a pair's rows.

Expected values are worked out by hand in each test, or read from the record itself."""
import numpy as np
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


class TestSilentChanging:
    def test_a_silent_and_a_changing_codon_are_counted_apart(self):
        a = "ATGCTAAAAGGGTTTCCCTAA"               # M L K G F P, stop
        b = "ATGCTGAAAGCGTTTCCCTAA"               # CTA>CTG silent (L); GGG>GCG changes G to A
        assert cybrid_dna.silent_changing(a, a) == (0, 0, 6)
        assert cybrid_dna.silent_changing(a, b) == (1, 1, 6)

    def test_the_mitochondrial_code_decides_what_is_silent(self):
        # ATA and ATG are both Met, TGA and TGG both Trp in the vertebrate mitochondrial code
        assert cybrid_dna.silent_changing("ATATGACCC", "ATGTGGCCC") == (2, 0, 3)

    def test_two_changed_bases_in_one_changing_codon_count_twice_and_a_gap_or_unknown_is_not_compared(self):
        assert cybrid_dna.silent_changing("ATGAAACCC", "ATGGGACCC") == (0, 2, 3)
        a = "ATGCTAAAAGGGTTTCCCTGGCATCAAGAA"
        assert cybrid_dna.silent_changing(a, a[:12] + a[15:])[2] == 9
        assert cybrid_dna.silent_changing(a, a[:12] + "NTT" + a[15:]) == (0, 0, 9)

    def test_codons_drop_a_final_stop_and_a_partial_end(self):
        assert cybrid_dna.codons("ATGAAATAA") == ["ATG", "AAA"]
        assert cybrid_dna.codons("ATGAAAT") == ["ATG", "AAA"]


class TestControlRegion:
    al = cybrid_dna.aligner()
    region = "ACGTTGCAAGGCTTAACCGGTTAACGTCAGTCAAGGTTCC" * 3

    def side(self, seq, blocks=None, start=100):
        return {"dloop": (seq, start), "blocks": blocks or {}}

    def test_not_assessed_without_an_annotation(self):
        got = cybrid_dna.control_region(self.side(self.region), {"dloop": None, "blocks": {}}, self.al)
        assert all(np.isnan(v) for v in got.values())

    def test_the_whole_region_and_the_blocks_of_the_host(self):
        donor = self.region[:10] + "T" + self.region[11:50] + "A" + self.region[51:]     # places 10 and 50 changed
        assert self.region[10] != "T" and self.region[50] != "A"
        blocks = {"CSB1": (145, 155), "CSB2": (160, 170), "CSB3": (180, 190)}            # places 45..54, 60..69, 80..89
        got = cybrid_dna.control_region(self.side(self.region, blocks), self.side(donor), self.al)
        assert got["nt_dloop_d"] == 2 and got["nt_dloop_f"] == 2 / len(self.region) and got["nt_csb_d"] == 1
        none = cybrid_dna.control_region(self.side(self.region), self.side(donor), self.al)
        assert none["nt_dloop_d"] == 2 and np.isnan(none["nt_csb_d"])

    def test_a_block_outside_the_region_is_refused(self):
        blocks = {"CSB1": (50, 60), "CSB2": (160, 170), "CSB3": (180, 190)}
        with pytest.raises(ValueError):
            cybrid_dna.control_region(self.side(self.region, blocks), self.side(self.region), self.al)


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

    def test_the_chinese_hamsters_two_trna_genes_are_read_on_the_minus_strand(self, features):
        # the record gives tRNA-Asn and tRNA-Cys on the plus strand (RECORD_FIXES): read as listed, each is as close
        # to the mouse's as the golden hamster's is, and far from it when reverse-complemented back
        al = cybrid_dna.aligner()
        mouse, chinese, golden = (features[cybrid.GENOMES[s]] for s in ("mus-musculus", "cricetulus-griseus", "mesocricetus-auratus"))
        for product in ("tRNA-Asn", "tRNA-Cys"):
            k = (product, 0)
            assert mouse["trna_strand"][k] == chinese["trna_strand"][k] == golden["trna_strand"][k] == -1
            ours = cybrid_dna.differences(mouse["trna"][k], chinese["trna"][k], al)[0]
            as_given = cybrid_dna.differences(mouse["trna"][k], str(Seq(chinese["trna"][k]).reverse_complement()), al)[0]
            near = cybrid_dna.differences(mouse["trna"][k], golden["trna"][k], al)[0]
            assert ours <= near + 5 and as_given >= 3 * ours

    def test_the_chinese_hamsters_control_region_is_the_stretch_after_trna_pro(self, features):
        seq, start = features[cybrid.GENOMES["cricetulus-griseus"]]["dloop"]
        assert (start, len(seq)) == (15403, 881)
        row = cybrid_dna.pair_rows(features, "mus-musculus", "cricetulus-griseus")[0]
        golden = cybrid_dna.pair_rows(features, "mus-musculus", "mesocricetus-auratus")[0]
        assert np.isfinite(row["nt_dloop_d"]) and abs(row["nt_dloop_d"] - golden["nt_dloop_d"]) < 0.2 * golden["nt_dloop_d"]

    def test_a_trna_gene_on_different_strands_in_the_two_records_is_refused(self, features):
        import copy
        f = copy.deepcopy(features)
        acc = cybrid.GENOMES["mus-spretus"]
        f[acc]["trna_strand"][("tRNA-Asn", 0)] = 1
        with pytest.raises(ValueError, match="different strands"):
            cybrid_dna.pair_rows(f, "mus-musculus", "mus-spretus")

    def test_a_species_against_itself_scores_nothing(self, features):
        for r in cybrid_dna.pair_rows(features, "mus-musculus", "mus-musculus"):
            assert all(r[s] == 0 for s in cybrid_dna.SCORES)

    def test_the_human_control_region_is_read_across_the_records_end_on_the_plus_strand(self, features):
        from Bio import SeqIO
        rec = next(r for r in SeqIO.parse(str(GB), "genbank") if r.id == "NC_012920.1")
        seq, start = features["NC_012920.1"]["dloop"]
        assert start == 16023 and seq == str(rec.seq[16023:] + rec.seq[:576]).upper()

    def test_the_mouse_blocks_lie_in_its_control_region_and_silent_and_changing_add_up(self, features):
        m = features[cybrid.GENOMES["mus-musculus"]]
        assert set(m["blocks"]) == set(cybrid_dna.CSB)
        rows = cybrid_dna.pair_rows(features, "mus-musculus", "mus-spretus")
        assert rows[0]["nt_csb_d"] <= rows[0]["nt_dloop_d"]
        for r in rows:
            assert r["nt_syn_d"] + r["nt_nonsyn_d"] <= r["nt_cds_d"] and r["nt_syn_d"] > r["nt_nonsyn_d"] > 0
        human = cybrid_dna.pair_rows(features, "homo-sapiens", "pan-troglodytes")
        assert np.isnan(human[0]["nt_csb_d"]) and human[0]["nt_dloop_d"] > 0

    def test_a_pair_is_the_same_both_ways_and_farther_donors_differ_more(self, features):
        ab = cybrid_dna.pair_rows(features, "homo-sapiens", "pan-troglodytes")
        ba = cybrid_dna.pair_rows(features, "pan-troglodytes", "homo-sapiens")
        far = cybrid_dna.pair_rows(features, "homo-sapiens", "pongo-abelii")
        assert [r["complex"] for r in ab] == list(cybrid.COMPLEX_GENES)
        for x, y, z in zip(ab, ba, far):
            assert x["nt_cds_d"] == y["nt_cds_d"] and x["nt_rrna_d"] == y["nt_rrna_d"]
            assert 0 < x["nt_cds_d"] < z["nt_cds_d"]
            assert 0 < x["nt_rrna_f"] < z["nt_rrna_f"] < 1
