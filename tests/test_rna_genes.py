"""H37 part A, step 1: the 22 tRNA and 2 rRNA genes of the records (`longevity/rna_genes.py`).

The records are real: the mitogenomes of the cybrid hosts and donors in `data/input/cybrid/records.gb` (`mitoage download
cybrid`), human, three apes, rat, mouse and a hamster. Every expected value is computed from the record's own text at run
time. The edge cases are the real record with one feature block removed, stripped of its names, or duplicated.
"""
import gzip
import hashlib
import re
import shutil

import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.build2 import genome as gm
from mitoage_pipelines.build2 import pins
from mitoage_pipelines.longevity import rna_genes as rg

RECORDS = paths.CYBRID_GB / "records.gb"
needs_records = pytest.mark.skipif(not RECORDS.is_file(), reason="needs `mitoage download cybrid`")
HUMAN = "NC_012920"
NAMES = ["NC_012920", "NC_001643", "NC_001644", "NC_011120", "NC_001665", "NC_005089", "NC_007936"]


@pytest.fixture(scope="module")
def records() -> dict[str, str]:
    texts = [t.lstrip() + "\n//\n" for t in RECORDS.read_text().split("\n//\n") if t.strip()]
    out = {re.search(r"^VERSION\s+(\S+)\.\d+", t, re.M).group(1): t for t in texts}
    return {k: out[k] for k in NAMES}


def blocks(text: str) -> tuple[str, list[str], str]:
    """A record split into the text before FEATURES, the feature blocks, and the text from ORIGIN on."""
    head, rest = text.split("\nFEATURES", 1)
    feat, tail = rest.split("\nORIGIN", 1)
    lines = ("FEATURES" + feat).split("\n")
    out, cur = [], []
    for line in lines[1:]:
        if re.match(r"^ {5}\S", line) and cur:
            out.append("\n".join(cur)); cur = []
        cur.append(line)
    out.append("\n".join(cur))
    return head + "\n" + lines[0], out, "\nORIGIN" + tail


def join(head: str, feats: list[str], tail: str) -> str:
    return head + "\n" + "\n".join(feats) + tail


def trnas(text: str) -> list[str]:
    return [b for b in blocks(text)[1] if re.match(r"^ {5}tRNA ", b)]


def anticodon(block: str) -> tuple[str, str]:
    m = re.search(r"aa:(\w+),seq:(\w+)", block)
    return m.group(1), m.group(2)


def expected_name(block: str) -> str:
    """The gene a tRNA block is, read from its anticodon, independently of build2/genome.py."""
    aa, seq = anticodon(block)
    iso = {("Leu", "taa"): "UUR", ("Leu", "tag"): "CUN", ("Ser", "tga"): "UCN", ("Ser", "gct"): "AGY"}
    return f"tRNA-{aa}({iso[(aa, seq)]})" if aa in ("Leu", "Ser") else f"tRNA-{aa}"


def feature_sequence(text: str, block: str) -> str:
    """The bases of a tRNA block in its own sense, from its location and the ORIGIN."""
    rec = gm.parse_genbank(text)
    loc = re.match(r"^ {5}\S+\s+(\S+)", block).group(1)
    comp = loc.startswith("complement")
    a, b = map(int, re.findall(r"\d+", loc))
    piece = rec.sequence[a - 1:b]
    return gm.reverse_complement(piece) if comp else piece


class TestExtract:
    @needs_records
    def test_every_trna_is_the_bases_of_its_location_in_its_own_sense(self, records):
        text = records[HUMAN]
        got = rg.extract_record(text, "Chordata")
        bl = trnas(text)
        assert len(bl) == 22
        assert sorted(expected_name(b) for b in bl) == sorted(rg.TRNA_GENES)
        for b in bl:
            assert got["sequences"][expected_name(b)] == feature_sequence(text, b)
        assert any(re.search(r"complement", b.split("\n")[0]) for b in bl)     # minus-strand genes are in the record
        assert set(got["sequences"]) == set(rg.RNA_GENES)
        assert got["unresolved"] == {} and set(got["copies"].values()) == {1}

    @needs_records
    def test_the_rrnas_are_the_two_rrna_features(self, records):
        text = records[HUMAN]
        rec = gm.parse_genbank(text)
        rr = [f for f in rec.features if f.key == "rRNA"]
        assert len(rr) == 2
        got = rg.extract_record(text, "Chordata")["sequences"]
        assert sorted([got["rRNA-12S"], got["rRNA-16S"]], key=len) == sorted([gm.extract(rec.sequence, f.segments).upper() for f in rr], key=len)
        assert len(got["rRNA-12S"]) < len(got["rRNA-16S"])

    @needs_records
    def test_a_record_missing_a_trna_leaves_that_gene_out_and_counts_it(self, records):
        text = records[HUMAN]
        head, feats, tail = blocks(text)
        gone = [b for b in feats if re.match(r"^ {5}tRNA ", b) and anticodon(b)[0] == "Trp"]
        assert len(gone) == 1
        cut = join(head, [b for b in feats if b is not gone[0]], tail)
        got = rg.extract_record(cut, "Chordata")
        assert "tRNA-Trp" not in got["sequences"] and len(got["sequences"]) == 23
        full = rg.extract_record(text, "Chordata")["sequences"]
        assert all(got["sequences"][g] == full[g] for g in got["sequences"])
        c = rg.census({"a": rg.extract_record(text, "Chordata"), "b": got})
        row = c.set_index("gene").loc["tRNA-Trp"]
        assert (row.n_gene, row.n_not_annotated, row.n_unresolved) == (1, 1, 0)

    @needs_records
    def test_leu_and_ser_without_names_are_placed_by_gene_order_in_a_chordate_and_unresolved_elsewhere(self, records):
        text = records[HUMAN]
        head, feats, tail = blocks(text)

        def strip(b: str) -> str:
            if re.match(r"^ {5}tRNA ", b) and anticodon(b)[0] in ("Leu", "Ser"):
                b = re.sub(r"\n {21}/(anticodon|codon_recognized)=[^\n]*", "", b)
                b = re.sub(r'/gene="[^"]*"', '/gene="TRNX"', b)
            return b
        bare = join(head, [strip(b) for b in feats], tail)
        full = rg.extract_record(text, "Chordata")["sequences"]
        chordate = rg.extract_record(bare, "Chordata")
        assert chordate["unresolved"] == {}
        for g in ("tRNA-Leu(UUR)", "tRNA-Leu(CUN)", "tRNA-Ser(UCN)", "tRNA-Ser(AGY)"):
            assert chordate["sequences"][g] == full[g]
        other = rg.extract_record(bare, "Arthropoda")
        assert other["unresolved"] == {"Leu": 2, "Ser": 2}
        assert not any(g in other["sequences"] for g in ("tRNA-Leu(UUR)", "tRNA-Leu(CUN)", "tRNA-Ser(UCN)", "tRNA-Ser(AGY)"))
        assert len(other["sequences"]) == 20
        c = rg.census({"human": rg.extract_record(text, "Chordata"), "bare": other}).set_index("gene")
        assert (c.loc["tRNA-Leu(UUR)", ["n_gene", "n_not_annotated", "n_unresolved", "n_unresolved_copies"]].tolist() == [1, 0, 1, 2])
        assert c.loc["tRNA-Ser(AGY)", "n_unresolved"] == 1 and c.loc["tRNA-Ala", "n_unresolved"] == 0

    @needs_records
    def test_a_duplicated_gene_keeps_the_first_copy_in_genome_order_and_is_counted(self, records):
        text = records[HUMAN]
        head, feats, tail = blocks(text)
        val = next(b for b in feats if re.match(r"^ {5}tRNA ", b) and anticodon(b)[0] == "Val")
        early = re.sub(r"^( {5}tRNA\s+)\S+", r"\g<1>1..70", val)          # a second copy placed before every other gene
        got = rg.extract_record(join(head, feats + [early], tail), "Chordata")
        assert got["copies"]["tRNA-Val"] == 2
        assert got["sequences"]["tRNA-Val"] == gm.parse_genbank(text).sequence[:70]
        late = re.sub(r"^( {5}tRNA\s+)\S+", r"\g<1>16000..16069", val)
        got2 = rg.extract_record(join(head, feats + [late], tail), "Chordata")
        assert got2["sequences"]["tRNA-Val"] == feature_sequence(text, val)
        c = rg.census({"a": got, "b": rg.extract_record(text, "Chordata")}).set_index("gene")
        assert c.loc["tRNA-Val", "n_duplicated"] == 1 and c.loc["tRNA-Val", "n_gene"] == 2

    @needs_records
    def test_a_trna_the_build_cannot_name_is_no_gene(self, tmp_path, records):
        text = records[HUMAN]
        head, feats, tail = blocks(text)
        trp = next(b for b in feats if re.match(r"^ {5}tRNA ", b) and anticodon(b)[0] == "Trp")
        odd = re.sub(r"aa:Trp", "aa:Sec", trp)
        odd = re.sub(r'/product="[^"]*"', '/product="tRNA-Sec"', odd)
        odd = re.sub(r'/gene="[^"]*"', '/gene="TRNU"', odd)
        got = rg.extract_record(join(head, [b for b in feats if b is not trp] + [odd], tail), "Chordata")
        assert "tRNA-Trp" not in got["sequences"] and len(got["sequences"]) == 23


class TestCensus:
    @needs_records
    def test_the_counts_of_a_gene_add_up_and_follow_the_records(self, records):
        per = {k: rg.extract_record(t, "Chordata") for k, t in records.items()}
        c = rg.census(per)
        assert list(c.gene) == rg.RNA_GENES
        assert (c.n_gene + c.n_not_annotated + c.n_unresolved == c.n_species).all()
        for g in rg.RNA_GENES:
            assert c.set_index("gene").loc[g, "n_gene"] == sum(g in r["sequences"] for r in per.values())
        assert (c.n_gene == len(records)).all()          # all seven records annotate all 24 genes


class TestFeatures:
    def test_gc_and_skews_are_the_build_statistics_of_the_bases_and_gaps_are_not_bases(self):
        s = "GGCCAATTTACGN"
        f = rg.gene_features("GG--CC-AATTT-ACGN")
        a, c, g, t = (s.count(x) for x in "ACGT")
        assert f["length"] == len(s) and f["n_other"] == 1
        assert f["gc_percent"] == round(100 * (g + c) / (a + c + g + t), 3)
        assert f["gc_skew"] == round((g - c) / (g + c), 5) and f["at_skew"] == round((a - t) / (a + t), 5)
        n = a + c + g + t
        assert (f["A_share"], f["C_share"], f["G_share"], f["T_share"]) == (a / n, c / n, g / n, t / n)

    def test_a_sequence_with_no_resolved_base_has_missing_shares_not_zeros(self):
        f = rg.gene_features("NNN---")
        assert f["length"] == 3 and f["A_share"] is None and f["gc_skew"] is None

    @needs_records
    def test_a_real_gene_agrees_with_the_build_statistics_of_its_region(self, records):
        text = records[HUMAN]
        g = gm.annotate(text, 2, True)
        feat = next(f for f in g.features if f.name == "tRNA-Leu(UUR)")
        st = gm.make_stats(0, "x", gm.count_bases(gm.extract(g.record.sequence, feat.segments)), "genbank")
        got = rg.gene_features(rg.extract_record(text, "Chordata")["sequences"]["tRNA-Leu(UUR)"])
        assert (got["gc_percent"], got["gc_skew"], got["at_skew"], got["length"]) == (st["gc_percent"], st["gc_skew"], st["at_skew"], st["length_bp"])


class TestInformativeColumns:
    @needs_records
    def test_counts_occupied_columns_with_two_bases_each_in_two_rows(self, records):
        real = [rg.extract_record(t, "Chordata")["sequences"]["tRNA-Val"] for t in list(records.values())[:6]]
        n = min(map(len, real))
        rows = {f"s{i}": s[:n] for i, s in enumerate(real)}
        rows["gapped"] = "-" * n                                            # a row of gaps changes no occupancy of the others
        rows["x"] = "N" * n
        # by hand, column by column, over the six real rows
        from collections import Counter
        want = 0
        for j in range(n):
            col = [r[j] for r in rows.values()]
            cnt = Counter(ch for ch in col if ch in "ACGT")
            want += (sum(cnt.values()) >= len(col) / 2) and sum(v >= 2 for v in cnt.values()) >= 2
        assert rg.informative_columns(rows) == want
        assert rg.informative_columns({k: v for k, v in rows.items() if k in {f"s{i}" for i in range(6)}}) >= want
        assert 0 < want < n

    def test_an_all_gap_column_a_singleton_and_a_sparse_column_are_not_informative(self):
        rows = {"a": "A-AC", "b": "A-AC", "c": "G-CT", "d": "G-C-", "e": "T-C-", "f": "A-GN"}
        # col 0: A,A,G,G,T,A : two bases with >= 2 rows -> informative; col 1: all gaps; col 2: A,A,C,C,C,G -> informative;
        # col 3: C,C,T only 3 of 6 rows have a base (exactly half) and C twice, T once -> one base with two rows only
        assert rg.informative_columns(rows) == 2

    def test_a_column_most_rows_lack_is_not_occupied_even_when_two_bases_are_each_seen_twice(self):
        rows = {"a": "A", "b": "A", "c": "C", "d": "C", **{f"g{i}": "-" for i in range(6)}}
        assert rg.informative_columns(rows) == 0                      # 4 of 10 rows have a base
        rows["e"], rows["f"] = "A", "C"
        assert rg.informative_columns(rows) == 1                      # 6 of 12: half the rows, the floor itself

    def test_fewer_than_four_rows_give_none_and_a_subset_counts_only_its_rows(self):
        rows = {"a": "AC", "b": "AC", "c": "GT", "d": "GT", "e": "AC"}
        assert rg.informative_columns({k: rows[k] for k in "abc"}) == 0
        assert rg.informative_columns(rows) == 2
        assert rg.informative_columns(rows, {"a", "b", "c", "d"}) == 2
        assert rg.informative_columns(rows, {"a", "b", "e", "c"}) == 0


@pytest.mark.skipif(shutil.which("mafft") is None, reason="needs MAFFT on PATH")
class TestAlign:
    @needs_records
    def test_each_row_is_its_input_in_slug_order_and_a_second_run_gives_the_same_text(self, records):
        seqs = {s: rg.extract_record(t, "Chordata")["sequences"]["tRNA-Val"] for s, t in records.items()}
        shuffled = {k: seqs[k] for k in sorted(seqs, reverse=True)}
        a, version = rg.align_genes({"tRNA-Val": shuffled, "tRNA-One": {"x": "ACGT", "y": "ACGA"}}, 2)
        b, _ = rg.align_genes({"tRNA-Val": seqs, "tRNA-One": {"x": "ACGT", "y": "ACGA"}}, 1)
        assert a == b
        assert a["tRNA-One"] == ""                               # fewer than three sequences: not aligned, said by an empty text
        rows = [r.split("\n") for r in a["tRNA-Val"].split(">")[1:]]
        assert [r[0] for r in rows] == sorted(seqs)
        assert len({len(r[1]) for r in rows}) == 1
        assert all(r[1].replace("-", "") == seqs[r[0]] for r in rows)
        assert version.startswith("v7.")


class TestLoadRecords:
    def _setup(self, root, records, in_release, monkeypatch, status=None):
        """A release of the records `in_release` (gzipped, its md5 made the pin) and a build naming all NAMES."""
        root.mkdir()
        rel = root / "mito.gz"
        with gzip.open(rel, "wt") as f:
            f.write("".join(records[n] for n in in_release))
        monkeypatch.setattr(pins, "REFSEQ_MD5", hashlib.md5(rel.read_bytes()).hexdigest())
        versions = [re.search(r"^VERSION\s+(\S+)", records[n], re.M).group(1) for n in NAMES]
        sp = pd.DataFrame({"slug": [f"sp-{n}" for n in NAMES], "refseq": versions, "phylum": "Chordata",
                           "record_status": status or "new"})
        build = root / "build"
        build.mkdir()
        sp.to_parquet(build / "species.parquet")
        return rel, build

    @needs_records
    def test_the_records_come_back_in_slug_order_with_their_text(self, tmp_path, records, monkeypatch):
        rel, build = self._setup(tmp_path / "a", records, NAMES[::-1], monkeypatch)
        got = rg.load_records(build, rel, tmp_path)
        assert [r["slug"] for r in got] == sorted(f"sp-{n}" for n in NAMES)
        assert all(r["text"] == records[r["slug"][3:]] for r in got)

    @needs_records
    def test_a_release_that_is_not_the_pinned_one_is_refused(self, tmp_path, records, monkeypatch):
        rel, build = self._setup(tmp_path / "a", records, NAMES, monkeypatch)
        monkeypatch.setattr(pins, "REFSEQ_MD5", "0" * 32)
        with pytest.raises(ValueError, match="not RefSeq release"):
            rg.load_records(build, rel, tmp_path)

    @needs_records
    def test_a_live_species_whose_record_is_not_in_the_release_is_refused_not_dropped(self, tmp_path, records, monkeypatch):
        rel, build = self._setup(tmp_path / "a", records, NAMES[:-1], monkeypatch)
        with pytest.raises(ValueError, match="1 live species have a record that is not in the release"):
            rg.load_records(build, rel, tmp_path)

    @needs_records
    def test_a_retired_species_is_left_out_even_when_its_record_is_missing(self, tmp_path, records, monkeypatch):
        rel, build = self._setup(tmp_path / "a", records, NAMES[:-1], monkeypatch, status=["new"] * (len(NAMES) - 1) + ["retired"])
        assert len(rg.load_records(build, rel, tmp_path)) == len(NAMES) - 1

    @needs_records
    def test_a_contig_only_record_takes_its_sequence_from_the_download_or_is_refused(self, tmp_path, records, monkeypatch):
        stub = records["NC_001665"].split("\nORIGIN")[0] + "\nCONTIG      join(NC_001665.2:1..16313)\n//\n"
        monkeypatch.setitem(records, "NC_001665", stub)
        rel, build = self._setup(tmp_path / "a", records, NAMES, monkeypatch)
        with pytest.raises(FileNotFoundError, match="CONTIG"):
            rg.load_records(build, rel, tmp_path / "nowhere")
        real = next(t.lstrip() + "\n//\n" for t in RECORDS.read_text().split("\n//\n") if re.search(r"^VERSION\s+NC_001665", t, re.M))
        c = tmp_path / "contig"
        c.mkdir()
        (c / f"{re.search(r'^VERSION\s+(\S+)', real, re.M).group(1)}.gb").write_text(real)
        got = {r["slug"]: r["text"] for r in rg.load_records(build, rel, c)}
        assert got["sp-NC_001665"] == real and got["sp-NC_005089"] != real


def test_the_22_trna_genes_are_the_20_amino_acids_with_two_isoacceptors_of_leu_and_ser():
    assert len(rg.TRNA_GENES) == 22 == len(set(rg.TRNA_GENES))
    assert sorted(g.split("(")[0] for g in rg.TRNA_GENES).count("tRNA-Leu") == 2
    assert {g.split("(")[0][5:] for g in rg.TRNA_GENES} == set(gm.AA3.values())
    assert rg.RNA_GENES == rg.TRNA_GENES + ["rRNA-12S", "rRNA-16S"]
    assert {rg.file_name(g) for g in rg.RNA_GENES}.__len__() == 24 and "(" not in "".join(map(rg.file_name, rg.RNA_GENES))
