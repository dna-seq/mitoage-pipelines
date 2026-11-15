"""Does aligning by structure change the columns a sequence aligner gives?

Per gene, over the species that have a predicted structure (`structure fold`): FoldMason builds the structural
multiple alignment, MAFFT L-INS-i the sequence one, on the same sequences. The two are compared as sets of
aligned residue pairs (sum-of-pairs agreement): for each pair of species, the residue pairs the structural
alignment puts in one column, and how many of them the sequence alignment also puts in one column; at most 500
species pairs per gene, drawn with seed 0. The same is reported for the conserved core alone: the columns of
the sequence alignment that at least 90 per cent of the species fill and at least 90 per cent of those carry
the modal residue. Both tools' run times are kept, as the cost.

Writes structural_msa.tsv. Ported from the cluster script `struct/smsa.py` of 2026-11-15 without change of
method.
"""
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .. import external
from ..seqio import read_fasta
from .fold import GENES


def residue_pairs(aln: dict[str, str], a: str, b: str) -> set[tuple[int, int]]:
    """The (0-based, ungapped) residue index pairs that the alignment puts in the same column."""
    ia = ib = -1; out = set()
    for x, y in zip(aln[a], aln[b]):
        if x != "-": ia += 1
        if y != "-": ib += 1
        if x != "-" and y != "-": out.add((ia, ib))
    return out


def core_columns(rows: list[str], min_fill: float = 0.9, min_modal: float = 0.9) -> list[int]:
    M = np.array([list(r) for r in rows]); core = []
    for j in range(M.shape[1]):
        col = M[:, j]; ok = col != "-"
        if ok.mean() < min_fill:
            continue
        _, n = np.unique(col[ok], return_counts=True)
        if n.max() / ok.sum() >= min_modal:
            core.append(j)
    return core


def foldmason_names(aln: dict[str, str]) -> dict[str, str]:
    """FoldMason names a row after its file (acinonyx-jubatus.pdb_A): back to the species slug."""
    return {k.replace(".pdb", "").split("_")[0]: v for k, v in aln.items()}


def main(pdbs: Path, build: Path, work: Path, out: Path, genes: list[str] | None = None, threads: int = 8) -> Path:
    mafft, foldmason = external.find("mafft", env="MAFFT"), external.find("foldmason")
    out.mkdir(parents=True, exist_ok=True); work.mkdir(parents=True, exist_ok=True)
    rows = []
    for g in genes or GENES:
        files = sorted((pdbs / g).glob("*.pdb"))
        if len(files) < 10:
            print(f"{g}: {len(files)} structures, skipped"); continue
        sp = [f.stem for f in files]
        seqs = read_fasta(build / "alignments" / f"{g}.aa.fasta.gz")
        fa = work / f"{g}.faa"
        fa.write_text("".join(f">{s}\n{seqs[s].replace('-', '')}\n" for s in sp))
        t0 = time.time()
        with open(work / f"{g}.seq.afa", "w") as f:
            external.run([mafft, "--localpair", "--maxiterate", "1000", "--anysymbol", "--quiet", "--thread", str(threads), fa], stdout=f)
        t_seq = time.time() - t0
        t0 = time.time()
        with open(work / f"{g}.fm.log", "w") as log:
            r = subprocess.run([foldmason, "easy-msa", *map(str, files), str(work / f"{g}.fm"), str(work / f"fm_{g}"),
                                "--report-mode", "0"], stdout=log, stderr=log)
        t_str = time.time() - t0
        fmf = next(iter(sorted(work.glob(f"{g}.fm*_aa.fa"))), None)
        if r.returncode != 0 or fmf is None:
            raise RuntimeError(f"{g}: FoldMason failed, see {work / f'{g}.fm.log'}")
        seq_aln, str_aln = read_fasta(work / f"{g}.seq.afa"), foldmason_names(read_fasta(fmf))
        sp = [s for s in sp if s in seq_aln and s in str_aln]
        core = core_columns([seq_aln[s] for s in sp])
        core_res = {}
        for s in sp:
            c = np.cumsum([ch != "-" for ch in seq_aln[s]]) - 1
            core_res[s] = {c[j] for j in core if seq_aln[s][j] != "-"}
        tot = agree = ctot = cagree = 0
        rng = np.random.default_rng(0)
        idx = [(a, b) for i, a in enumerate(sp) for b in sp[i + 1:]]
        if len(idx) > 500:
            idx = [idx[i] for i in rng.choice(len(idx), 500, replace=False)]
        for a, b in idx:
            ps, pq = residue_pairs(str_aln, a, b), residue_pairs(seq_aln, a, b)
            tot += len(ps); agree += len(ps & pq)
            cs = {(x, y) for x, y in ps if x in core_res[a] and y in core_res[b]}
            ctot += len(cs); cagree += len(cs & pq)
        rows.append(dict(gene=g, species=len(sp), columns_seq=len(seq_aln[sp[0]]), core_columns=len(core),
                         sp_agreement=agree / max(tot, 1), core_agreement=cagree / max(ctot, 1),
                         pairs_scored=len(idx), seconds_mafft=t_seq, seconds_foldmason=t_str))
        print(f"{g}: {len(sp)} species, agreement {agree / max(tot, 1):.4f} overall, {cagree / max(ctot, 1):.4f} in the "
              f"core ({len(core)} columns); mafft {t_seq:.0f} s, foldmason {t_str:.0f} s", flush=True)
    res = pd.DataFrame(rows)
    table = out / "structural_msa.tsv"
    res.to_csv(table, sep="\t", index=False, float_format="%.4f")
    if not res.empty:
        print(f"\nmedian agreement over genes: {res.sp_agreement.median():.4f} overall, "
              f"{res.core_agreement.median():.4f} in the conserved core")
    print("wrote", table)
    return table
