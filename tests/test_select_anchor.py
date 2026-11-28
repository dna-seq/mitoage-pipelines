"""The human anchor of `mitoport/select.py`: which human protein every species' isoform is chosen nearest to.

The rows below are the human RefSeq proteins of two real NCBI mammalian ortholog sets downloaded on 2026-11-28
(`mitoage download orthologs --gene CARS2 --gene SARS2`), in NCBI's own file order and with their own lengths. CARS2
is the case the option exists for: its first protein in file order is a 302-residue isoform and its canonical one is
564, so the two anchors disagree by a factor of nearly two. SARS2 is the ordinary case, two proteins 2 residues apart,
where the choice cannot matter much.

`main` itself is not covered here: it reads two FASTA files and a build parquet, and the choice this module is about
is the one line these tests hold.
"""
import pytest

from mitoage_pipelines.mitoport import select

# (accession, organism, gene id, sequence); only the lengths matter, so the sequences are that many residues
CARS2 = [("NP_001339181.1", 302), ("NP_001339182.1", 351), ("NP_078813.1", 564), ("XP_006720016.1", 451),
         ("XP_047286565.1", 302), ("XP_054230933.1", 302)]
SARS2 = [("NP_001139373.1", 520), ("NP_060297.1", 518)]


def prots(rows, organism="Homo sapiens", gene_id="79587"):
    return [(a, organism, gene_id, "M" * n) for a, n in rows]


def test_first_is_ncbi_file_order():
    assert select.human_length(prots(CARS2), "first") == 302
    assert select.human_length(prots(SARS2), "first") == 520


def test_longest_is_the_longest_refseq_protein():
    assert select.human_length(prots(CARS2), "longest") == 564
    assert select.human_length(prots(SARS2), "longest") == 520


def test_the_two_anchors_differ_on_the_case_the_option_is_for():
    """If this ever stops differing, the option has no case and the tests above are vacuous."""
    assert select.human_length(prots(CARS2), "first") != select.human_length(prots(CARS2), "longest")


def test_only_refseq_proteins_of_the_human_count():
    """An XP_ model protein is never the anchor, and another species' protein is not either, however long."""
    rows = prots([("XP_999999.1", 5000)]) + prots(CARS2)
    assert select.human_length(rows, "longest") == 564
    mouse = [("NP_111111.1", "Mus musculus", "99", "M" * 9000)]
    assert select.human_length(prots(CARS2) + mouse, "longest") == 564


def test_a_set_with_no_human_refseq_protein_is_refused():
    with pytest.raises(ValueError, match="no human RefSeq protein"):
        select.human_length(prots([("XP_006720016.1", 451)]), "longest")
    with pytest.raises(ValueError, match="no human RefSeq protein"):
        select.human_length([("NP_111111.1", "Mus musculus", "99", "M" * 400)], "first")


def test_an_unknown_anchor_is_refused_rather_than_defaulted():
    with pytest.raises(ValueError, match="expected 'first' or 'longest'"):
        select.human_length(prots(CARS2), "canonical")


def test_the_default_is_the_old_behaviour():
    """Every run before 2026-11-28 anchored on the first human protein; the default must still do that, or a
    published number moves without anyone asking for it."""
    import inspect
    assert inspect.signature(select.main).parameters["anchor"].default == "first"
