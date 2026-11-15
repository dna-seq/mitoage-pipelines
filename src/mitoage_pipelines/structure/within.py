"""All-against-all structural superposition of the proteins of the human respiratory complexes.

Chains come from solved human structures (mmCIF, from `download pdb`), tried in the order given; the first
structure that carries a gene supplies its chain. A chain is assigned to a mitochondrially encoded protein by
sequence (identity of at least 0.8 to MitoAge's human protein, length within 0.5 to 1.6 of it), because the
mmCIF descriptions name the nuclear accessory subunits of complex I in the same words as the encoded ones;
the nuclear subunits of complex IV (5Z62 only) are named by description. Every pair is superimposed with
US-align; the TM-score normalised by the longer chain is the one to read (normalising by a short chain
inflates it), and the one normalised by the shorter chain, the aligned length and the RMSD are kept.

The mitochondrially encoded set is not one family: ND2, ND4 and ND5 descend from an Mrp-type antiporter and
should superimpose on each other, the rest should not. That is the control of this step.

Writes within_species.tsv and chains.tsv, and one PDB file per gene's chain. Ported from the cluster script
`struct/within.py` of 2026-11-14; the published run tried 5XTD 5Z62 8H9S 5XTE 5XTH 5XTI 6ZQN 7AJD. Unlike that
script it refuses to run without US-align instead of falling back to RMSD-only superposition.
"""
import itertools
import re
import subprocess
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBIO, MMCIFParser, Select
from Bio.PDB.Polypeptide import protein_letters_3to1

from .. import external
from ..genes import MT_GENES
from ..mitoport.burial import gene_of as nuclear_gene_of
from ..mitoport.contacts import entity_descriptions
from ..seqio import read_fasta

CANDIDATES = "5XTD 5Z62 8H9S 5XTE 5XTH 5XTI 6ZQN 7AJD".split()
MT13 = "ND1 ND2 ND3 ND4 ND4L ND5 ND6 CYTB COX1 COX2 COX3 ATP6 ATP8".split()     # the published row order
assert sorted(MT13) == sorted(MT_GENES)


class One(Select):
    def __init__(self, cid): self.cid = cid
    def accept_chain(self, ch): return ch.id == self.cid
    def accept_residue(self, r): return r.id[0] == " "


def parse_usalign(stdout: str) -> dict | None:
    """TM-scores (shorter- and longer-chain normalised), aligned length and RMSD from US-align's report."""
    tm = [float(m) for m in re.findall(r"TM-score=\s*([0-9.]+)", stdout)]
    al = re.search(r"Aligned length=\s*(\d+), RMSD=\s*([0-9.]+)", stdout)
    if len(tm) < 2 or not al:
        return None
    return dict(tm_short=max(tm[:2]), tm_long=min(tm[:2]), aligned=int(al.group(1)), rmsd=float(al.group(2)))


def main(build: Path, cif_dir: Path, candidates: list[str], chains_dir: Path, out: Path) -> Path:
    usalign = external.find("USalign")
    out.mkdir(parents=True, exist_ok=True); chains_dir.mkdir(parents=True, exist_ok=True)
    human = {g: read_fasta(build / "alignments" / f"{g}.aa.fasta.gz")["homo-sapiens"].replace("-", "") for g in MT13}
    aligner = PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1,
                              substitution_matrix=substitution_matrices.load("BLOSUM62"))

    def mt_gene(seq: str):
        best, ident = None, 0.0
        for g, h in human.items():
            if not (0.5 < len(seq) / len(h) < 1.6):
                continue
            a = aligner.align(seq, h)[0]
            n = sum(sum(seq[s0 + k] == h[h0 + k] for k in range(s1 - s0)) for (s0, s1), (h0, h1) in zip(*a.aligned))
            f = n / min(len(seq), len(h))
            if f > ident:
                best, ident = g, f
        return (best, ident) if ident >= 0.8 else (None, ident)

    io = PDBIO(); have = {}
    for pdb in candidates:
        cif = cif_dir / f"{pdb}.cif"
        if not cif.is_file():
            raise FileNotFoundError(f"{cif} is missing: run `mitoage download pdb {' '.join(candidates)}`")
        desc_of = entity_descriptions(cif)
        st = MMCIFParser(QUIET=True).get_structure(pdb, str(cif))[0]
        found = []
        for ch in st:
            desc = desc_of.get(ch.id, "?")
            res = [r for r in ch if r.id[0] == " "]
            if len(res) < 30:
                continue
            seq = "".join(protein_letters_3to1.get(r.get_resname().upper(), "X") for r in res)
            g, ident = mt_gene(seq)
            if g is None:
                g, ident = (nuclear_gene_of(desc), float("nan")) if pdb == "5Z62" else (None, float("nan"))
            if not g or g in have:
                continue
            io.set_structure(st); io.save(str(chains_dir / f"{g}.pdb"), select=One(ch.id))
            have[g] = dict(gene=g, pdb=pdb, chain=ch.id, residues=len(res), mt=int(g in MT13),
                           identity_to_human=round(ident, 3) if ident == ident else "", description=desc)
            found.append(f"{g}({ch.id},{len(res)}" + (f",id {ident:.2f})" if ident == ident else ")"))
        print(f"{pdb}: {' '.join(found) if found else 'nothing new'}", flush=True)
    pd.DataFrame(have.values()).sort_values("gene").to_csv(out / "chains.tsv", sep="\t", index=False)
    print(f"\nmitochondrially encoded proteins found: {sorted(g for g in have if g in MT13)}")
    print(f"missing: {[g for g in MT13 if g not in have]}")

    rows = []
    for a, b in itertools.combinations(sorted(have), 2):
        r = subprocess.run([usalign, str(chains_dir / f"{a}.pdb"), str(chains_dir / f"{b}.pdb")],
                           capture_output=True, text=True, check=True)
        res = parse_usalign(r.stdout)
        if res is None:
            raise RuntimeError(f"US-align gave no TM-score for {a} against {b}:\n{r.stdout[-500:]}")
        rows.append(dict(a=a, b=b, **res))
    res = pd.DataFrame(rows)
    table = out / "within_species.tsv"
    res.to_csv(table, sep="\t", index=False, float_format="%.4f")
    mt = res[res.a.isin(MT13) & res.b.isin(MT13)]
    print(f"\nmitochondrially encoded pairs: {len(mt)}; above TM-score 0.5 normalised by the longer chain:")
    print(mt[mt.tm_long > 0.5].sort_values("tm_long", ascending=False).to_string(index=False))
    print("wrote", table, out / "chains.tsv")
    return table
