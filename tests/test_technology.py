"""H8's control: reading a record's sequencing technology and sorting it into a class."""
import polars as pl
import pytest

from mitoage_pipelines.longevity import motif_pgls, technology

RECORD = """LOCUS       {acc}              16000 bp    DNA     circular VRT 01-JAN-2024
VERSION     {acc}
COMMENT     {comment}
FEATURES             Location/Qualifiers
     source          1..16000
//
"""

ASSEMBLY = """##Assembly-Data-START##
            Assembly Method       :: MitoHiFi v. 2.2
            Sequencing Technology :: {tech}
            ##Assembly-Data-END##"""


def record(acc, comment):
    return RECORD.format(acc=acc, comment=comment)


class TestStatedTechnology:
    def test_reads_the_value_of_the_structured_comment(self):
        assert technology.stated_technology(record("OX1.1", ASSEMBLY.format(tech="PacBio Sequel"))) == "PacBio Sequel"

    def test_joins_a_value_wrapped_over_two_lines(self):
        wrapped = ASSEMBLY.format(tech="Illumina HiSeq 2500; Oxford\n            Nanopore MinION")
        assert technology.stated_technology(record("OX1.1", wrapped)) == "Illumina HiSeq 2500; Oxford Nanopore MinION"

    def test_stops_at_the_next_field(self):
        c = "##Assembly-Data-START##\n            Sequencing Technology :: Sanger dideoxy sequencing\n            Coverage              :: 30x\n            ##Assembly-Data-END##"
        assert technology.stated_technology(record("OX1.1", c)) == "Sanger dideoxy sequencing"

    def test_a_record_without_the_field_states_nothing(self):
        assert technology.stated_technology(record("NC_1.1", "REVIEWED REFSEQ: derived from AB000001.")) is None

    def test_a_record_without_a_comment_states_nothing(self):
        assert technology.stated_technology("LOCUS x\nFEATURES             Location/Qualifiers\n//\n") is None

    def test_an_empty_value_states_nothing(self):
        assert technology.stated_technology(record("OX1.1", ASSEMBLY.format(tech=""))) is None


class TestClassify:
    @pytest.mark.parametrize("text,cls", [
        ("PacBio Sequel", "long_read"), ("Oxford Nanopore PromethION", "long_read"), ("ONT", "long_read"),
        ("Illumina NovaSeq; PacBio HiFi", "long_read"),            # a hybrid assembly was spanned by long reads
        ("Illumina HiSeq 2000", "short_read"), ("454", "short_read"), ("IonTorrent", "short_read"), ("BGISEQ-500", "short_read"),
        ("Sanger dideoxy sequencing", "sanger"), ("ABI 3730", "sanger"),
        ("", "unstated"), (None, "unstated"),
    ])
    def test_classes(self, text, cls):
        assert technology.classify(text) == cls

    def test_ont_inside_a_word_is_not_nanopore(self):
        # "ont" as a substring would make these long-read
        assert technology.classify("Illumina; contig assembly, front end") == "short_read"
        assert technology.classify("Montage of amplicons") == "unstated"

    def test_a_text_naming_no_platform_is_unstated_never_short_read(self):
        assert technology.classify("next generation sequencing") == "unstated"


class TestTable:
    def build(self, tmp_path):
        data = tmp_path / "build"; data.mkdir()
        pl.DataFrame({"slug": ["a", "b", "c", "d"], "refseq": ["NC_1.1", "NC_2.1", "NC_3.1", None]}).write_parquet(data / "species.parquet")
        folder = tmp_path / "tech"; (folder / "refseq").mkdir(parents=True); (folder / "source").mkdir()
        # a: the source states long reads; b: only the RefSeq record states it; c: RefSeq record present, source named but not fetched
        (folder / "refseq" / "NC_1.1.gb").write_text(record("NC_1.1", "PROVISIONAL REFSEQ: This record is derived from OX1."))
        (folder / "source" / "OX1.2.gb").write_text(record("OX1.2", ASSEMBLY.format(tech="PacBio Sequel II")))
        (folder / "refseq" / "NC_2.1.gb").write_text(record("NC_2.1", "REVIEWED REFSEQ: identical to AB2.\n            " + ASSEMBLY.format(tech="Illumina MiSeq")))
        (folder / "source" / "AB2.1.gb").write_text(record("AB2.1", "no assembly data"))
        (folder / "refseq" / "NC_3.1.gb").write_text(record("NC_3.1", "PROVISIONAL REFSEQ: derived from KX3.1."))
        return data, folder

    def test_source_first_then_refseq_then_unstated(self, tmp_path):
        data, folder = self.build(tmp_path)
        t = {r["slug"]: r for r in technology.table(data, folder).iter_rows(named=True)}
        assert (t["a"]["technology_class"], t["a"]["stated_in"], t["a"]["source"]) == ("long_read", "source", "OX1")
        assert (t["b"]["technology_class"], t["b"]["stated_in"]) == ("short_read", "refseq")
        assert (t["c"]["technology_class"], t["c"]["source"], t["c"]["source_record"]) == ("unstated", "KX3.1", False)
        assert (t["d"]["technology_class"], t["d"]["refseq_record"]) == ("unstated", False)
        assert set(t) == {"a", "b", "c", "d"}                      # no species is dropped

    def test_main_writes_the_table_the_motif_test_reads(self, tmp_path):
        data, folder = self.build(tmp_path)
        out = tmp_path / "out"
        technology.main(data, folder, out)
        assert motif_pgls.long_read_slugs(None, out / "technology.tsv") == {"a"}

    def test_main_refuses_without_the_download(self, tmp_path):
        data, _ = self.build(tmp_path)
        with pytest.raises(SystemExit):
            technology.main(data, tmp_path / "nothing", tmp_path / "out")

    def test_where_both_records_state_one_the_source_is_kept(self, tmp_path):
        data, folder = self.build(tmp_path)
        (folder / "refseq" / "NC_1.1.gb").write_text(record("NC_1.1", "PROVISIONAL REFSEQ: derived from OX1.\n            " + ASSEMBLY.format(tech="Illumina MiSeq")))
        a = technology.table(data, folder).filter(pl.col("slug") == "a").row(0, named=True)
        assert (a["technology"], a["technology_class"], a["stated_in"]) == ("PacBio Sequel II", "long_read", "source")
