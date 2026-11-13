"""Fourfold-site counts on a tiny codon alignment; the expected counts are derived here from the
vertebrate mitochondrial code, not from the list of codon families the code under test carries."""
import gzip
from itertools import product

from Bio.Data import CodonTable

from mitoage_pipelines.genes import MT_HEAVY_STRAND
from mitoage_pipelines.longevity import fourfold

ALN = {
    #        CTA GTC TCG CCT ACA GCC CGG GGT | not fourfold: ATA TTT AAA | gap, N third base, lower case
    "a-one": "CTAGTCTCGCCTACAGCCCGGGGT" "ATATTTAAA" "---" "CTN" "gcc",
    "b-two": "CTTCTTCTTCTTCTC---CTG------ATGATGATG" "TGA" "GTA" "ccc",
    "c-none": "ATGAAATTTTAA------",
}


def expected_counts(seq: str) -> dict[str, int]:
    """Third bases of the codons whose first two bases code one amino acid whatever the third (table 2)."""
    fwd = CodonTable.unambiguous_dna_by_id[2].forward_table
    families = {a + b for a, b in product("ACGT", repeat=2) if len({fwd.get(a + b + c) for c in "ACGT"}) == 1
                and fwd.get(a + b + "A") is not None}
    cnt = dict.fromkeys("ACGT", 0)
    s = seq.upper()
    for i in range(0, len(s) - 2, 3):
        if s[i:i + 2] in families and s[i + 2] in "ACGT":
            cnt[s[i + 2]] += 1
    return cnt


def test_count_fourfold_matches_the_genetic_code():
    got = fourfold.count_fourfold(ALN)
    want = {k: expected_counts(v) for k, v in ALN.items()}
    assert got == want
    assert sum(want["a-one"].values()) == 9          # 8 families + the lower-case GCC; CTN is not counted
    assert want["c-none"] == dict.fromkeys("ACGT", 0)


def test_counts_accumulate_over_genes():
    once = fourfold.count_fourfold(ALN)
    twice = fourfold.count_fourfold(ALN, fourfold.count_fourfold(ALN))
    assert twice == {k: {b: 2 * n for b, n in v.items()} for k, v in once.items()}


def test_main_writes_shares_of_the_summed_counts(tmp_path):
    data = tmp_path / "build"; (data / "alignments").mkdir(parents=True)
    for i, g in enumerate(MT_HEAVY_STRAND):          # a different slice of the species in each gene
        rows = {k: v for j, (k, v) in enumerate(ALN.items()) if (i + j) % 3 or k == "a-one"}
        with gzip.open(data / "alignments" / f"{g}.codon.fasta.gz", "wt") as f:
            f.write("".join(f">{k} some description\n{v[:21]}\n{v[21:]}\n" for k, v in rows.items()))
    fourfold.main(data, tmp_path / "out")
    lines = (tmp_path / "out" / "fourfold.tsv").read_text().splitlines()
    assert lines[0].split("\t") == ["slug", "A4", "C4", "T4", "fourfold_sites"]
    got = {l.split("\t")[0]: l.split("\t")[1:] for l in lines[1:]}
    assert list(got) == sorted(ALN)
    for j, (k, v) in enumerate(ALN.items()):
        genes_with = sum(1 for i in range(len(MT_HEAVY_STRAND)) if (i + j) % 3 or k == "a-one")
        c = {b: n * genes_with for b, n in expected_counts(v).items()}
        n = max(1, sum(c.values()))
        assert got[k] == [f"{c['A'] / n:.5f}", f"{c['C'] / n:.5f}", f"{c['T'] / n:.5f}", str(n)]
