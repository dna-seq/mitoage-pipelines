"""Gene lists, the download record, the aligner-agreement scores and the branch table, each on a small real input."""
import json

import numpy as np

from mitoage_pipelines import download, genes
from mitoage_pipelines.longevity import compare
from mitoage_pipelines.mitoport import erc


def test_gene_lists_are_disjoint_and_complete():
    lists = [genes.COMPLEX_IV_NUCLEAR, genes.COMPLEX_III_NUCLEAR, genes.CYTOCHROME_C, genes.COMPLEX_II,
             genes.MATRIX, genes.CYTOSOLIC_CONTROLS]
    assert genes.NUCLEAR == [g for l in lists for g in l]
    assert len(set(genes.NUCLEAR)) == len(genes.NUCLEAR) == sum(len(l) for l in lists)
    assert {g for v in genes.GROUPS.values() for g in v} == set(genes.NUCLEAR) - set(genes.CYTOSOLIC_CONTROLS)
    assert set(genes.MT_GENES) - set(genes.MT_HEAVY_STRAND) == {"ND6"}
    assert {m.removeprefix("MT-") for ms in genes.MT_OF_COMPLEX.values() for m in ms} <= set(genes.MT_GENES)


def test_manifest_appends_with_counts_from_the_folder(tmp_path):
    folder = tmp_path / "input" / "thing"; (folder / "sub").mkdir(parents=True)
    contents = {"a.txt": b"abc", "sub/b.bin": b"\x00" * 10}
    for name, blob in contents.items(): (folder / name).write_bytes(blob)
    (folder / ".cache").mkdir(); (folder / ".cache" / "lock").write_bytes(b"not counted")
    manifest = tmp_path / "input" / "MANIFEST.json"
    first = download.append_manifest(manifest, "test", "one", folder)
    (folder / "c.txt").write_bytes(b"12345")
    second = download.append_manifest(manifest, "test", "two", folder, revision="r1")
    records = json.loads(manifest.read_text())
    assert records == [first, second]
    assert (first["files"], first["bytes"]) == (len(contents), sum(map(len, contents.values())))
    assert (second["files"], second["bytes"]) == (first["files"] + 1, first["bytes"] + 5)
    assert second["revision"] == "r1" and first["downloaded_at"].endswith("+00:00")


def brute_force_pairs(a: dict[str, str], b: dict[str, str], ids: list[str]):
    """Per column of a: residue pairs it aligns, and how many of them b also puts in one column."""
    def col_of(aln):          # (row, k-th residue) -> column
        return {(i, k): c for i in ids for k, c in enumerate(c for c, ch in enumerate(aln[i]) if ch != "-")}
    ca, cb = col_of(a), col_of(b)
    ncol = len(a[ids[0]]); total, kept = [0] * ncol, [0] * ncol
    cells = sorted(ca)
    for n, p in enumerate(cells):
        for q in cells[n + 1:]:
            if p[0] != q[0] and ca[p] == ca[q]:
                total[ca[p]] += 1; kept[ca[p]] += cb[p] == cb[q]
    return total, kept


def test_aligner_agreement_counts_pairs():
    ids = ["x", "y", "z"]
    a = {"x": "AC-G", "y": "ACTG", "z": "A-TG"}       # the reference
    b = {"x": "ACG-", "y": "ACTG", "z": "AT-G"}       # x's G and z's T moved one column
    ia, ib = compare.residue_index(a, ids), compare.residue_index(b, ids)
    total, kept, exact = compare.score(ia, compare.column_of_residue(ib))
    want_total, want_kept = brute_force_pairs(a, b, ids)
    assert (total.tolist(), kept.tolist()) == (want_total, want_kept)
    assert want_total == [3, 1, 1, 3] and want_kept == [3, 1, 0, 1]   # col 3: the two T split; col 4: only y and z stay together
    assert exact.tolist() == [True, True, False, False]
    same_total, same_kept, same_exact = compare.score(ia, compare.column_of_residue(ia))
    assert np.array_equal(same_total, same_kept) and same_exact.all()


def test_branches_are_keyed_by_the_smaller_side(tmp_path):
    p = tmp_path / "t.treefile"
    p.write_text("((a:1,b:2):0.5,c:3,(d:4,e:5):0.25);")
    got = erc.branches(p)
    f = frozenset
    assert got == {f("a"): 1.0, f("b"): 2.0, f("c"): 3.0, f("d"): 4.0, f("e"): 5.0, f("ab"): 0.5, f("de"): 0.25}
    # the same unrooted tree written from another root gives the same table; a root split in two adds up
    p.write_text("((a:1,b:2):0.2,((d:4,e:5):0.25,c:3):0.3);")
    assert erc.branches(p) == got
