"""How buried each mitochondrial residue of complex IV is, and a non-contact control set matched on burial.

Solvent accessibility by Shrake-Rupley twice on the solved human enzyme: once in the whole complex, once
for each mitochondrial chain alone. Relative accessibility divides by the theoretical maxima of Tien et al.
2013. A residue is a contact when any heavy atom lies within 4.5 A of a nuclear-encoded subunit
(contacts.tsv, from `mitoport contacts`). Each contact residue is matched, within its own gene and without
replacement, to the nearest non-contact residue by accessibility in the isolated chain, the measure that does
not depend on the partner being present. Every residue is mapped to its column of the MitoAge alignment,
with that column's occupancy and identity to human across the species of the analysis: those with an
alignment, a tip on the dated tree and a lifespan record of usable quality.

Writes sites.tsv (one row per residue: group contact, control or other), chains.tsv and species.txt.
Ported from the cluster script `mitoport/iface/burial.py` of 2026-11-14 without change of method.
"""
import copy
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

from ..seqio import read_fasta
from .contacts import entity_descriptions

MTCHAIN = {"A": "COX1", "B": "COX2", "C": "COX3"}       # the chain ids of PDB 5Z62
MAXASA = dict(ALA=129, ARG=274, ASN=195, ASP=193, CYS=167, GLU=223, GLN=225, GLY=104, HIS=224,
              ILE=197, LEU=201, LYS=236, MET=224, PHE=240, PRO=159, SER=155, THR=172, TRP=285,
              TYR=263, VAL=174)
NUCNAME = [(r"subunit 4 isoform 1", "COX4I1"), (r"subunit 5A", "COX5A"), (r"subunit 5B", "COX5B"),
           (r"subunit 6A1", "COX6A1"), (r"subunit 6B1", "COX6B1"), (r"subunit 6C", "COX6C"),
           (r"subunit 7A2", "COX7A2"), (r"subunit 7B", "COX7B"), (r"subunit 7C", "COX7C"),
           (r"subunit 8A", "COX8A"), (r"NDUFA4", "NDUFA4")]


def gene_of(description: str) -> str | None:
    """The nuclear gene of complex IV an entity description names, else None."""
    for pat, g in NUCNAME:
        if re.search(pat, description, re.I):
            return g
    return None


def column_of_position(row: str) -> dict[int, int]:
    """{1-based residue position in the ungapped row: 0-based column of the alignment}."""
    out, k = {}, 0
    for c, ch in enumerate(row):
        if ch != "-":
            k += 1; out[k] = c
    return out


def tree_tips(newick: str) -> set[str]:
    return set(re.findall(r"[(,]([A-Za-z0-9\-\._]+):", newick))


def match_controls(group: pd.DataFrame) -> tuple[list, list]:
    """(contact index, control index) for one gene: each contact residue, in order of increasing isolated-chain
    accessibility, takes the nearest unused non-contact residue that has an alignment column. Ties in the
    distance go to the first in table order (pandas idxmin). When the pool runs out the remaining contacts go unmatched."""
    con = group[group.n_contacts > 0].sort_values("rasa_alone")     # numpy quicksort, as published: deterministic
    pool = group[(group.n_contacts == 0) & (group.aln_col > 0)].copy()
    taken = []
    for _, row in con.iterrows():
        if pool.empty:
            break
        j = (pool.rasa_alone - row.rasa_alone).abs().idxmin()
        taken.append(j); pool = pool.drop(j)
    return list(con.index), taken


def main(cif: Path, contacts: Path, build: Path, tree: Path, quality_ok: Path, out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    chain_desc = entity_descriptions(cif)
    model = MMCIFParser(QUIET=True).get_structure("x", str(cif))[0]
    for ch in model:                                     # waters, haems, lipids and ions are not protein
        for r in [r for r in ch if r.id[0] != " "]:
            ch.detach_child(r.id)
    desc = {ch.id: chain_desc.get(ch.id, "?") for ch in model}
    pd.DataFrame([{"chain": c, "description": v, "gene": MTCHAIN.get(c) or gene_of(v),
                   "is_mt": int(c in MTCHAIN), "residues": len(model[c])} for c, v in desc.items()]
                 ).to_csv(out / "chains.tsv", sep="\t", index=False)
    sr = ShrakeRupley()
    sr.compute(model, level="R")
    sasa_complex = {(ch.id, r.id): r.sasa for ch in model for r in ch}
    sasa_alone = {}
    for cid in MTCHAIN:
        m = copy.deepcopy(model)
        for ch in [ch for ch in m if ch.id != cid]:
            m.detach_child(ch.id)
        sr.compute(m, level="R")
        for r in m[cid]:
            sasa_alone[(cid, r.id)] = r.sasa

    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    recs, alns = [], {}
    for cid, gene in MTCHAIN.items():
        aln = read_fasta(build / "alignments" / f"{gene}.aa.fasta.gz"); alns[gene] = aln
        hrow = aln["homo-sapiens"]; hseq = hrow.replace("-", "")
        col_of_pos = column_of_position(hrow)
        res = list(model[cid])
        seq = "".join(protein_letters_3to1.get(r.get_resname().upper(), "X") for r in res)
        a = al.align(seq, hseq)[0]
        m = {s0 + k: h0 + k + 1 for (s0, s1), (h0, h1) in zip(*a.aligned) for k in range(s1 - s0)}
        print(f"{gene} chain {cid}: {len(seq)} residues, {len(m)} mapped, identity to the human protein "
              f"{sum(seq[i] == hseq[j - 1] for i, j in m.items()) / len(seq):.3f}")
        for i, r in enumerate(res):
            hp = m.get(i); mx = MAXASA.get(r.get_resname().upper(), np.nan)
            recs.append(dict(gene=gene, chain=cid, pdb_resnum=r.id[1], resname=r.get_resname(),
                             human_pos=hp or 0, human_aa=hseq[hp - 1] if hp else "",
                             aln_col=col_of_pos.get(hp, 0) + 1 if hp else 0,
                             rasa_complex=sasa_complex[(cid, r.id)] / mx, rasa_alone=sasa_alone[(cid, r.id)] / mx))
    sites = pd.DataFrame(recs)
    sites["d_rasa"] = sites.rasa_alone - sites.rasa_complex

    ct = pd.read_csv(contacts, sep="\t")
    nuc = ct[ct.partner_is_mt == 0].copy()
    nuc["partner_gene"] = [gene_of(x) or x for x in nuc.partner]
    per = nuc.groupby(["mt_gene", "pdb_resnum"]).agg(
        n_contacts=("partner_resnum", "size"),
        partners=("partner_gene", lambda s: ",".join(sorted(set(s))))).reset_index()
    sites = sites.merge(per.rename(columns={"mt_gene": "gene"}), on=["gene", "pdb_resnum"], how="left")
    sites["n_contacts"] = sites.n_contacts.fillna(0).astype(int); sites["partners"] = sites.partners.fillna("")

    species = sorted(set(alns["COX1"]) & tree_tips(tree.read_text()) & set(quality_ok.read_text().split()))
    (out / "species.txt").write_text("\n".join(species) + "\n")
    print(f"species with an alignment, a tip on the dated tree and a lifespan of usable quality: {len(species)}")
    occ, ident = [], []
    for g, grp in sites.groupby("gene"):
        M = np.array([list(alns[g][s]) for s in species])
        cols = grp.aln_col.to_numpy() - 1
        gap = (M[:, cols] == "-").mean(axis=0)
        same = np.array([(M[:, c] == alns[g]["homo-sapiens"][c]).mean() for c in cols])
        occ.append(pd.Series(1 - gap, index=grp.index)); ident.append(pd.Series(same, index=grp.index))
    sites["occupancy"] = pd.concat(occ); sites["identity_to_human"] = pd.concat(ident)

    sites["group"] = "other"
    for _, grp in sites.groupby("gene"):
        con, taken = match_controls(grp)
        sites.loc[con, "group"] = "contact"
        sites.loc[taken, "group"] = "control"
    table = out / "sites.tsv"
    sites.to_csv(table, sep="\t", index=False, float_format="%.4f")

    print("\nper gene and group: residues, median relative accessibility (isolated chain / in the complex), "
          "median column occupancy, median identity to human")
    for (g, grp), dd in sites.groupby(["gene", "group"]):
        print(f"  {g} {grp:8s} n={len(dd):3d}  rASA alone {dd.rasa_alone.median():.3f}  in complex "
              f"{dd.rasa_complex.median():.3f}  occupancy {dd.occupancy.median():.3f}  identity {dd.identity_to_human.median():.3f}")
    print("wrote", table, out / "chains.tsv", out / "species.txt")
    return table
