"""H37 part C's gene set (`mitoport/synthetases.py`): the synthetases charging the mitochondrial tRNAs, as NCBI gives them.

The fixture `data/ncbi_trna_ligase_genes.xml` is nine real NCBI gene records, trimmed to their track id, their Gene-ref
and their GeneOntology commentary verbatim, out of the efetch answer for the module's own search term on 2026-11-28.
They are the nine cases the rule has to tell apart, and each expectation below is read off the record's own text by a
path other than the code's: the description is matched with a plain string method, the GO terms are counted from the
raw bytes, and nothing is compared against a list of symbols kept in a test.

No test here touches the network. `search` and `fetch` are not covered: they are two URL calls, and what they return is
the fixture.
"""
import re
from pathlib import Path

import pytest

from mitoage_pipelines.mitoport import synthetases as syn

FIXTURE = Path(__file__).parent / "data" / "ncbi_trna_ligase_genes.xml"


@pytest.fixture(scope="module")
def rows():
    return syn.parse(FIXTURE.read_bytes())


@pytest.fixture(scope="module")
def by_symbol(rows):
    return {r["symbol"]: r for r in rows}


def test_every_record_is_read_once(rows):
    """One row per record, and nothing read twice: the count comes from the raw bytes, not from the parser."""
    assert len(rows) == FIXTURE.read_bytes().count(b"<Entrezgene>")
    assert len({r["symbol"] for r in rows}) == len(rows)
    assert all(r["gene_id"].isdigit() and r["description"] for r in rows)


def test_the_description_decides_mitochondrial(rows):
    """Every record whose NCBI description ends in ", mitochondrial" is `mitochondrial`, and no other is."""
    said = {r["symbol"] for r in rows if r["description"].endswith(", mitochondrial")}
    assert {r["symbol"] for r in rows if r["status"] == "mitochondrial"} == said
    assert said, "the fixture must hold at least one described-mitochondrial gene"


def test_a_mitochondrial_component_without_a_mitochondrial_process_is_ambiguous(rows):
    """The three-way rule, checked against the terms of each record rather than against a list of genes."""
    for r in rows:
        if r["status"] == "mitochondrial":
            continue
        components, processes = r["components"].split("; "), r["processes"].split("; ")
        has_component = any(t in syn.MITO_COMPONENTS for t in components)
        has_process = any(p.startswith("mitochondrial ") and "aminoacylation" in p for p in processes)
        expected = "other" if not has_component else ("dual" if has_process else "ambiguous")
        assert r["status"] == expected, f"{r['symbol']}: {r['status']} for components {components} processes {processes}"


def test_all_three_statuses_and_a_used_and_an_unused_gene_are_in_the_fixture(rows):
    """The cases the rule has to tell apart are all present, so the tests above are not vacuous."""
    assert {r["status"] for r in rows} >= {"mitochondrial", "dual", "ambiguous", "other"}


def test_the_amino_acid_comes_from_the_description_not_from_the_processes(by_symbol):
    """DARS2 is the record that makes this matter: NCBI's process term for it names asparagine and its description
    names aspartate. The code must follow the description."""
    dars2 = by_symbol["DARS2"]
    assert "asparaginyl" in dars2["processes"], "the fixture no longer holds the disagreeing record this test is for"
    assert dars2["description"].startswith("aspartyl-")
    assert (dars2["amino_acid"], dars2["code"]) == ("aspartate", "Asp")


def test_the_amino_acid_matches_the_prefix_of_every_description(rows):
    """Read independently: the prefix up to the first hyphen or space, looked up in the naming table."""
    for r in rows:
        prefix = re.split(r"[- ]", r["description"])[0].lower()
        expected = syn.ACYL.get(prefix) if r["description"].lower().startswith(prefix + "-tRNA synthetase".lower()) or \
            r["description"].lower().startswith(prefix + " trna synthetase") else None
        if expected is None:
            continue
        assert (r["amino_acid"], r["code"]) == expected, r["symbol"]


def test_a_description_with_two_amino_acids_resolves_to_none(by_symbol):
    """EPRS1 is `glutamyl-prolyl-tRNA synthetase 1`: two amino acids, so the code names neither."""
    eprs1 = by_symbol["EPRS1"]
    assert eprs1["description"].count("yl-") >= 2
    assert syn.amino_acid(eprs1["description"]) is None
    assert eprs1["code"] == ""


def test_a_description_without_a_hyphen_is_read(by_symbol):
    """WARS2's description is `tryptophanyl tRNA synthetase 2, mitochondrial`, with a space where the others have a
    hyphen; a parser that required the hyphen would drop a used gene."""
    wars2 = by_symbol["WARS2"]
    assert "tryptophanyl tRNA" in wars2["description"]
    assert (wars2["code"], wars2["status"]) == ("Trp", "mitochondrial")


def test_only_mitochondrial_and_dual_genes_are_used(rows):
    used, rest = syn.resolve(rows)
    assert {r["status"] for r in used.values()} <= set(syn.USED)
    assert all(r["status"] not in syn.USED for r in rest)
    assert len(used) + len(rest) == len(rows)
    assert {r["code"] for r in used.values()} == set(used)


def test_two_genes_for_one_amino_acid_are_refused(rows):
    """A set this code cannot resolve stops the run; it does not pick one."""
    first = next(r for r in rows if r["status"] in syn.USED)
    twin = dict(first, symbol=first["symbol"] + "X", gene_id="999999")
    with pytest.raises(ValueError, match=f"two genes for {first['code']}"):
        syn.resolve(rows + [twin])


def test_a_used_gene_with_no_amino_acid_is_refused(rows):
    broken = dict(next(r for r in rows if r["status"] in syn.USED), amino_acid="", code="")
    with pytest.raises(ValueError, match="names no amino acid"):
        syn.resolve([broken])


def test_every_amino_acid_without_a_gene_is_named_with_its_reason(rows):
    """The 20 amino acids are covered: each is either used or listed as not assessed, with a reason that names the
    candidate genes when there were any."""
    used, rest = syn.resolve(rows)
    missing = syn.not_assessed(used, rest)
    codes = {c for _, (_, c) in syn.ACYL.items()}
    assert set(used) | {c for c, _, _ in missing} == codes
    assert len(used) + len(missing) == len(codes)
    for code, name, why in missing:
        candidates = {r["symbol"] for r in rest if r["code"] == code}
        for s in candidates:
            assert s in why, f"{code}: {s} is a candidate and the reason does not name it"
        if not candidates:
            assert "no human gene" in why


def test_the_table_written_holds_a_row_for_every_amino_acid(tmp_path, monkeypatch):
    """The run end to end on the fixture's bytes, with no network: `main` reads the file it is given."""
    out = tmp_path / "synthetases.tsv"
    monkeypatch.setattr(syn, "search", lambda: pytest.fail("the run must not search when an answer is given"))
    monkeypatch.setattr(syn, "fetch", lambda ids: pytest.fail("the run must not fetch when an answer is given"))
    syn.main(out, FIXTURE)
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == syn.COLUMNS
    body = [dict(zip(syn.COLUMNS, line.split("\t"))) for line in lines[1:]]
    assert len(body) == len(syn.ACYL)
    assert {r["code"] for r in body} == {c for _, c in syn.ACYL.values()}
    for r in body:
        assert (r["status"] == "not assessed") == (r["symbol"] == "")
