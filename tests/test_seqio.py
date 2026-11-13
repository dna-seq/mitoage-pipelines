"""The FASTA reader on files the test writes: gzip and plain, several extensions, wrapped lines."""
import gzip

import pytest

from mitoage_pipelines.seqio import read_fasta, read_records, write_fasta

RECORDS = [  # name, description, sequence as wrapped lines
    ("homo-sapiens", "NP_004246.2 [organism=Homo sapiens] [GeneID=9377]", ["MLGAA", "lrrc-", "-SV"]),
    ("mus-musculus", "", ["--ACGT"]),
    ("zebra_fish.1", "x", ["A"]),
]
TEXT = "".join(f">{n}{' ' + d if d else ''}\n" + "\n".join(lines) + "\n" for n, d, lines in RECORDS)
EXPECTED = {n: "".join(lines) for n, _, lines in RECORDS}


@pytest.mark.parametrize("name", ["x.fasta.gz", "x.aa.fasta.gz", "x.fasta", "x.fa", "x.faa", "x.fna", "x.afa"])
def test_reads_gzip_and_plain(tmp_path, name):
    p = tmp_path / name
    if name.endswith(".gz"):
        with gzip.open(p, "wt") as f: f.write(TEXT)
    else:
        p.write_text(TEXT)
    got = read_fasta(p)
    assert got == EXPECTED                               # wrapped lines joined, case and gaps kept
    assert list(got) == [n for n, _, _ in RECORDS]       # file order
    assert read_records(p) == [(n, d, "".join(lines)) for n, d, lines in RECORDS]


def test_round_trip(tmp_path):
    p = tmp_path / "out.faa"
    write_fasta(p, EXPECTED)
    assert read_fasta(p) == EXPECTED
    assert list(read_fasta(p)) == list(EXPECTED)


def test_many_records_keep_their_order(tmp_path):
    seqs = {f"s{i:04d}": "ACDEFGHIKL"[i % 10:] + "-" * (i % 7) for i in range(3000)}
    p = tmp_path / "many.fasta.gz"
    with gzip.open(p, "wt") as f:
        f.write("".join(f">{k} n={i}\n{v}\n" for i, (k, v) in enumerate(seqs.items())))
    got = read_fasta(p)
    assert list(got.items()) == list(seqs.items())


def test_duplicate_name_is_refused(tmp_path):
    p = tmp_path / "dup.fasta"
    p.write_text(">a one\nAC\n>b\nGG\n>a two\nTT\n")
    with pytest.raises(ValueError, match="duplicate sequence name a"):
        read_fasta(p)
    assert [r[0] for r in read_records(p)] == ["a", "b", "a"]      # the record view keeps both


def test_empty_and_missing(tmp_path):
    p = tmp_path / "empty.fasta"
    p.write_text("")
    assert read_fasta(p) == {}
    with pytest.raises(FileNotFoundError):
        read_fasta(tmp_path / "absent.fasta.gz")
