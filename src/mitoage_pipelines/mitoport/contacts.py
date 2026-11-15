"""Contacts between the mitochondrially and the nuclear encoded subunits of a solved complex IV.

A residue pair is in contact when any two of their heavy atoms are within CUTOFF angstrom (4.5). Chains are
named by the entity description of the mmCIF file. Each mitochondrial residue is mapped to its position in
the human protein of the MitoAge alignment by a global pairwise alignment of the chain's sequence to it.

Writes contacts.tsv: one row per (mitochondrial residue, partner residue) pair in contact, partners of
either origin, with `partner_is_mt` saying which. The published run used the human enzyme PDB 5Z62.
Ported from the cluster script `mitoport/structure/contacts.py` of 2026-11-14 without change of method.
"""
import collections
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import protein_letters_3to1

from ..seqio import read_fasta

CUTOFF = 4.5
MTKEY = {"COX1": ("subunit 1", "subunit i "), "COX2": ("subunit 2", "subunit ii "), "COX3": ("subunit 3", "subunit iii")}
HEADER = ("mt_gene\tchain\tpdb_resnum\tpdb_res\thuman_pos\thuman_aa\tpartner_chain\tpartner\tpartner_is_mt"
          "\tpartner_resnum\tpartner_res\n")


def mt_gene(description: str) -> str | None:
    """COX1, COX2 or COX3 when an entity description names one of them, else None."""
    t = description.lower() + " "
    for g, keys in MTKEY.items():
        if any(k in t for k in keys) and "nadh" not in t:
            return g
    return None


def human_protein(build: Path, gene: str) -> str:
    return read_fasta(build / "alignments" / f"{gene}.aa.fasta.gz")["homo-sapiens"].replace("-", "")


def entity_descriptions(cif: Path) -> dict[str, str]:
    """{chain id: entity description} from the mmCIF dictionary."""
    d = MMCIF2Dict(str(cif))
    ent_desc = dict(zip(d["_entity.id"], d["_entity.pdbx_description"]))
    out = {}
    for ent, chains in zip(d["_entity_poly.entity_id"], d["_entity_poly.pdbx_strand_id"]):
        for c in chains.split(","):
            out[c.strip()] = ent_desc.get(ent, "?")
    return out


def main(cif: Path, build: Path, out: Path) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    chain_desc = entity_descriptions(cif)
    model = MMCIFParser(QUIET=True).get_structure(cif.stem, str(cif))[0]
    desc = {c.id: chain_desc.get(c.id, "?") for c in model}
    for c, v in desc.items():
        print(f"  chain {c}: {v} ({sum(1 for r in model[c] if r.id[0] == ' ')} residues)")
    mt_chain = {c: mt_gene(v) for c, v in desc.items() if mt_gene(v)}
    print("mitochondrial chains:", mt_chain)

    al = PairwiseAligner(); al.mode = "global"; al.open_gap_score = -10; al.extend_gap_score = -0.5
    pos = {}
    for c, g in mt_chain.items():
        res = [r for r in model[c] if r.id[0] == " "]
        seq = "".join(protein_letters_3to1.get(r.get_resname().upper(), "X") for r in res)
        h = human_protein(build, g)
        a = al.align(seq, h)[0]; m = {}
        for (s0, s1), (h0, h1) in zip(*a.aligned):
            for k in range(s1 - s0):
                m[s0 + k] = h0 + k + 1
        ident = sum(seq[i] == h[j - 1] for i, j in m.items()) / len(seq)
        print(f"  {g} chain {c}: {len(seq)} residues, identity to the human protein {ident:.3f}")
        pos[c] = {res[i].id: (j, h[j - 1]) for i, j in m.items()}

    atoms = [a for c in model for r in c if r.id[0] == " " for a in r if a.element != "H"]
    pairs = set()
    for a, b in NeighborSearch(atoms).search_all(CUTOFF):
        ra, rb = a.get_parent(), b.get_parent(); ca, cb = ra.get_parent().id, rb.get_parent().id
        if ca == cb:
            continue
        for (r1, c1), (r2, c2) in (((ra, ca), (rb, cb)), ((rb, cb), (ra, ca))):
            if c1 in mt_chain:
                pairs.add((c1, r1.id, r1.get_resname(), c2, r2.id[1], r2.get_resname()))

    table = out / "contacts.tsv"
    with open(table, "w") as f:
        f.write(HEADER)
        for c1, rid, rn, c2, r2, rn2 in sorted(pairs, key=lambda x: (mt_chain[x[0]], x[1][1], x[3], x[4], x)):   # x: a stable tie-break
            hp = pos[c1].get(rid, ("", ""))
            f.write(f"{mt_chain[c1]}\t{c1}\t{rid[1]}\t{rn}\t{hp[0]}\t{hp[1]}\t{c2}\t{desc[c2]}\t{int(c2 in mt_chain)}"
                    f"\t{r2}\t{rn2}\n")
    touching = collections.defaultdict(set)
    for c1, rid, rn, c2, r2, rn2 in pairs:
        if c2 not in mt_chain:
            touching[mt_chain[c1]].add(rid)
    for g in sorted(set(mt_chain.values())):
        tot = sum(1 for c in mt_chain if mt_chain[c] == g for r in model[c] if r.id[0] == " ")
        print(f"{g}: {len(touching[g])} of {tot} residues touch a nuclear subunit")
    print("wrote", table)
    return table
