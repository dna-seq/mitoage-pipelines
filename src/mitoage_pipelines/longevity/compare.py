"""Agreement between alignments of the same sequences made by different tools.

Takes alignment files named <GENE>.<tool>.afa; for each gene writes agreement.tsv (one row per gene and
ordered pair of tools) and columns.<GENE>.tsv (per column of the reference tool, the share of its
residue pairs each other tool also aligns).

sum-of-pairs(A, B): of the residue pairs A puts in one column, the share B also puts in one column.
column score(A, B): of A's columns with two residues or more, the share B reproduces exactly.
"""
from pathlib import Path

import numpy as np

from ..seqio import read_fasta

REF = "mafft_fast"


def residue_index(aln: dict[str, str], ids: list[str]) -> np.ndarray:
    """rows x columns: the running index of the residue in its sequence, -1 at a gap."""
    lengths = {len(aln[i]) for i in ids}
    if len(lengths) != 1:
        raise SystemExit(f"rows of different lengths: {sorted(lengths)[:5]}")
    m = np.array([np.frombuffer(aln[i].encode(), dtype=np.uint8) for i in ids])
    res = m != ord("-")
    idx = np.cumsum(res, axis=1) - 1
    idx[~res] = -1
    return idx


def column_of_residue(idx: np.ndarray) -> list[np.ndarray]:
    """per row: column of the k-th residue."""
    return [np.nonzero(row >= 0)[0] for row in idx]


def pairs(n: np.ndarray) -> np.ndarray:
    return n * (n - 1) // 2


def score(a: np.ndarray, b_cols: list[np.ndarray]):
    """Per column of A: (pairs A aligns, pairs B also aligns, column reproduced exactly)."""
    rows, cols = a.shape
    # where B puts each residue of A's cell
    in_b = np.full(a.shape, -1, dtype=np.int64)
    for r in range(rows):
        cell = a[r] >= 0
        in_b[r, cell] = b_cols[r][a[r, cell]]
    total = np.zeros(cols, dtype=np.int64)
    kept = np.zeros(cols, dtype=np.int64)
    exact = np.zeros(cols, dtype=bool)
    for c in range(cols):
        v = in_b[:, c]
        v = v[v >= 0]
        if v.size < 2:
            continue
        total[c] = pairs(np.int64(v.size))
        _, counts = np.unique(v, return_counts=True)
        kept[c] = pairs(counts.astype(np.int64)).sum()
        exact[c] = counts.size == 1
    return total, kept, exact


def by_gene(files: list[Path]) -> dict[str, dict[str, Path]]:
    """gene -> tool -> file, from names <GENE>.<tool>.afa; genes in order of first appearance, tools by file name."""
    out: dict[str, dict[str, Path]] = {}
    for p in files:
        parts = p.name.split(".")
        if len(parts) != 3 or parts[2] != "afa":
            raise SystemExit(f"{p}: an alignment is named <GENE>.<tool>.afa")
        out.setdefault(parts[0], {})
    for gene in out:
        for p in sorted((p for p in files if p.name.split(".")[0] == gene), key=lambda p: p.name):
            tool = p.name.split(".")[1]
            if tool in out[gene]:
                raise SystemExit(f"{gene}: two files for the tool {tool}")
            out[gene][tool] = p
    return out


def main(files: list[Path], out: Path, ref: str = REF) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rows_out = ["gene\ta\tb\tsequences\tcolumns_a\tcolumns_b\tsum_of_pairs\tcolumn_score\tsum_of_pairs_occupied\tcolumn_score_occupied"]
    for gene, tools_files in by_gene(files).items():
        alns = {t: {k: v.upper() for k, v in read_fasta(p).items()} for t, p in tools_files.items()}
        if ref not in alns:
            raise SystemExit(f"{gene}: no {ref} alignment")
        ids = sorted(alns[ref])
        for t, a in alns.items():
            if sorted(a) != ids:
                raise SystemExit(f"{gene}: {t} has other sequences than {ref}")
            for i in ids:
                if a[i].replace("-", "") != alns[ref][i].replace("-", ""):
                    raise SystemExit(f"{gene}: {t} changed the residues of {i}")
        idx = {t: residue_index(a, ids) for t, a in alns.items()}
        cols = {t: column_of_residue(m) for t, m in idx.items()}
        percol = {}
        for ta in idx:
            for tb in idx:
                if ta == tb:
                    continue
                total, kept, exact = score(idx[ta], cols[tb])
                occupied = (idx[ta] >= 0).sum(axis=0) * 2 >= len(ids)
                scored = total > 0
                so = scored & occupied
                rows_out.append("\t".join(map(str, [
                    gene, ta, tb, len(ids), idx[ta].shape[1], idx[tb].shape[1],
                    f"{kept.sum() / total.sum():.4f}", f"{exact[scored].mean():.4f}",
                    f"{kept[so].sum() / total[so].sum():.4f}", f"{exact[so].mean():.4f}",
                ])))
                if ta == ref:
                    percol[tb] = (total, kept, occupied)
        tools = sorted(percol)
        with open(out / f"columns.{gene}.tsv", "w") as f:
            f.write("column\toccupied\t" + "\t".join(tools) + "\n")
            n = idx[ref].shape[1]
            for c in range(n):
                vals = [f"{percol[t][1][c] / percol[t][0][c]:.3f}" if percol[t][0][c] else "" for t in tools]
                f.write(f"{c + 1}\t{int(percol[tools[0]][2][c])}\t" + "\t".join(vals) + "\n")
        print(gene, "done", list(alns), flush=True)
    (out / "agreement.tsv").write_text("\n".join(rows_out) + "\n")
