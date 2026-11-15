"""The solved human complexes the contact maps are read from, and which chain of each file is which gene.

Every protein entity of an entry is named here by its UniProt entry (the mmCIF's `_struct_ref`, database UNP) and
mapped to the gene symbol its ortholog set is downloaded under. A file holding an entity this table does not name,
or an entity without a UniProt reference, is refused rather than guessed: the chain ids of a deposited entry are
the depositors' choice. Each mitochondrially encoded gene has one reference chain, the one whose residues are mapped
and scored; where an entry holds two copies (complex III is a dimer), the other copy is a mitochondrial partner and
its residues are never counted as nuclear contacts. Checked by hand against the four entries on 2026-11-15
(`_entity.pdbx_description`, `_struct_ref.db_code`, `_entity_poly.pdbx_strand_id`):

- 5Z62, human complex IV: COX1-3 on chains A, B, C; eleven nuclear subunits, one chain each.
- 5XTD, human complex I: ND1 s, ND2 i, ND3 j, ND4 r, ND4L k, ND5 l, ND6 m; 37 nuclear entities (the acyl carrier
  protein NDUFAB1 twice, chains G and X).
- 5XTE, human complex III, a dimer: cytochrome b on J and V (J is the reference); nine nuclear genes in two copies,
  the Rieske protein UQCRFS1 as two entities (its cleaved presequence, B and O, and the mature protein, C and P).
- 8H9S, human ATP synthase: ATP6 on N, ATP8 on Q; the c-ring of eight ATP5MC1 chains (1 to 8); alpha and beta in
  three copies; the inhibitor ATP5IF1 (J) is bound in this state and is mapped like the rest.
"""
from dataclasses import dataclass, field
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from .. import paths
from ..genes import COMPLEX_I_NUCLEAR, COMPLEX_III_NUCLEAR, COMPLEX_IV_NUCLEAR, COMPLEX_V_NUCLEAR, MT_GENES


@dataclass(frozen=True)
class Complex:
    name: str                       # I, III, IV, V
    pdb: str                        # the entry read
    mt: dict[str, str]              # mitochondrial gene -> its reference chain, in the order the gene is processed
    uniprot: dict[str, str]         # UniProt entry name -> gene, every protein entity of the entry
    nuclear: list[str]              # the nuclear genes of the entry, in the order their rows are written
    anchor: str                     # nuclear gene whose species make H9's "with nuclear partners" set (`mitoport iface`)
    groups: list[tuple[str, list[str]]] = field(default_factory=list)   # other groups of `mitoport iface`, in order

    @property
    def cif(self) -> Path:
        return paths.PDB / f"{self.pdb}.cif"

    @property
    def work(self) -> Path:
        """contacts.tsv, sites.tsv, control_pairs.tsv, chains.tsv, species.txt; complex IV's keeps its first place."""
        return paths.MITOPORT_IFACE if self.name == "IV" else paths.MITOPORT / f"iface-{self.name}"

    @property
    def out(self) -> Path:
        """coupling.tsv and lifespan.tsv of `mitoport iface`; complex IV's keeps its first place."""
        return paths.IFACE_OUT if self.name == "IV" else paths.MITOPORT_OUT / f"complex-{self.name}"

    @property
    def partners_out(self) -> Path:
        """partners.tsv and partners_species.tsv of `mitoport partners`; complex IV's keeps its first place."""
        return paths.MITOPORT_OUT if self.name == "IV" else paths.MITOPORT_OUT / f"complex-{self.name}"

    @property
    def partner_group(self) -> str:
        return f"complex {self.name} partner"


# --- the entries ---------------------------------------------------------------------------------
_IV = {"COX1_HUMAN": "COX1", "COX2_HUMAN": "COX2", "COX3_HUMAN": "COX3", "COX41_HUMAN": "COX4I1",
       "COX5A_HUMAN": "COX5A", "COX5B_HUMAN": "COX5B", "CX6A1_HUMAN": "COX6A1", "CX6B1_HUMAN": "COX6B1",
       "COX6C_HUMAN": "COX6C", "CX7A2_HUMAN": "COX7A2", "COX7B_HUMAN": "COX7B", "COX7C_HUMAN": "COX7C",
       "COX8A_HUMAN": "COX8A", "NDUA4_HUMAN": "NDUFA4"}
_I = {"H9PGF0_HUMAN": "ND1", "Q4GRX1_HUMAN": "ND2", "B9EE38_HUMAN": "ND3", "B2XJG5_HUMAN": "ND4",
      "V9JN72_HUMAN": "ND4L", "X5BVZ3_HUMAN": "ND5", "Q8HAX7_HUMAN": "ND6",
      "NDUV1_HUMAN": "NDUFV1", "NDUV2_HUMAN": "NDUFV2", "NDUV3_HUMAN": "NDUFV3",
      "NDUS1_HUMAN": "NDUFS1", "NDUS2_HUMAN": "NDUFS2", "NDUS3_HUMAN": "NDUFS3", "NDUS4_HUMAN": "NDUFS4",
      "NDUS5_HUMAN": "NDUFS5", "NDUS6_HUMAN": "NDUFS6", "NDUS7_HUMAN": "NDUFS7", "NDUS8_HUMAN": "NDUFS8",
      "NDUA1_HUMAN": "NDUFA1", "NDUA2_HUMAN": "NDUFA2", "NDUA3_HUMAN": "NDUFA3", "NDUA5_HUMAN": "NDUFA5",
      "NDUA6_HUMAN": "NDUFA6", "NDUA7_HUMAN": "NDUFA7", "NDUA8_HUMAN": "NDUFA8", "NDUA9_HUMAN": "NDUFA9",
      "NDUAA_HUMAN": "NDUFA10", "NDUAB_HUMAN": "NDUFA11", "NDUAC_HUMAN": "NDUFA12", "NDUAD_HUMAN": "NDUFA13",
      "ACPM_HUMAN": "NDUFAB1",
      "NDUB1_HUMAN": "NDUFB1", "NDUB2_HUMAN": "NDUFB2", "NDUB3_HUMAN": "NDUFB3", "NDUB4_HUMAN": "NDUFB4",
      "NDUB5_HUMAN": "NDUFB5", "NDUB6_HUMAN": "NDUFB6", "NDUB7_HUMAN": "NDUFB7", "NDUB8_HUMAN": "NDUFB8",
      "NDUB9_HUMAN": "NDUFB9", "NDUBA_HUMAN": "NDUFB10", "NDUBB_HUMAN": "NDUFB11",
      "NDUC1_HUMAN": "NDUFC1", "NDUC2_HUMAN": "NDUFC2"}
_III = {"CYB_HUMAN": "CYTB", "QCR1_HUMAN": "UQCRC1", "QCR2_HUMAN": "UQCRC2", "CY1_HUMAN": "CYC1",
        "UCRI_HUMAN": "UQCRFS1", "QCR7_HUMAN": "UQCRB", "QCR8_HUMAN": "UQCRQ", "QCR6_HUMAN": "UQCRH",
        "QCR9_HUMAN": "UQCR10", "QCR10_HUMAN": "UQCR11"}
_V = {"ATP6_HUMAN": "ATP6", "ATP8_HUMAN": "ATP8", "ATPA_HUMAN": "ATP5F1A", "ATPB_HUMAN": "ATP5F1B",
      "ATPG_HUMAN": "ATP5F1C", "ATPD_HUMAN": "ATP5F1D", "ATP5E_HUMAN": "ATP5F1E", "AT5F1_HUMAN": "ATP5PB",
      "ATP5H_HUMAN": "ATP5PD", "ATP5J_HUMAN": "ATP5PF", "ATPO_HUMAN": "ATP5PO", "AT5G1_HUMAN": "ATP5MC1",
      "ATP5I_HUMAN": "ATP5ME", "ATPK_HUMAN": "ATP5MF", "ATP5L_HUMAN": "ATP5MG", "ATP68_HUMAN": "ATP5MJ",
      "ATIF1_HUMAN": "ATP5IF1"}

# the other protein groups of `mitoport iface`, as H5 ran them for complex IV (whose own nuclear genes are the partner)
_H5_GROUPS = [("complex IV nuclear", COMPLEX_IV_NUCLEAR), ("complex III nuclear", COMPLEX_III_NUCLEAR)]


def _groups(own: list[str]) -> list[tuple[str, list[str]]]:
    return [(n, g) for n, g in _H5_GROUPS if g != own]


COMPLEXES = {
    "I": Complex("I", "5XTD", {"ND1": "s", "ND2": "i", "ND3": "j", "ND4": "r", "ND4L": "k", "ND5": "l", "ND6": "m"},
                 _I, COMPLEX_I_NUCLEAR, "NDUFS2", _groups(COMPLEX_I_NUCLEAR)),
    "III": Complex("III", "5XTE", {"CYTB": "J"}, _III, COMPLEX_III_NUCLEAR, "UQCRC1", _groups(COMPLEX_III_NUCLEAR)),
    "IV": Complex("IV", "5Z62", {"COX1": "A", "COX2": "B", "COX3": "C"}, _IV, COMPLEX_IV_NUCLEAR, "COX4I1",
                  _groups(COMPLEX_IV_NUCLEAR)),
    "V": Complex("V", "8H9S", {"ATP6": "N", "ATP8": "Q"}, _V, COMPLEX_V_NUCLEAR, "ATP5F1A", _groups(COMPLEX_V_NUCLEAR)),
}


def get(name: str) -> Complex:
    if name not in COMPLEXES:
        raise SystemExit(f"unknown complex {name!r}: one of {', '.join(COMPLEXES)}")
    return COMPLEXES[name]


def chain_genes(cif: Path, cx: Complex) -> dict[str, str]:
    """{chain id: gene} for every protein chain of `cif`, read from its UniProt references and checked against `cx`:
    refused when an entity has no UniProt reference or one `cx` does not name, when a gene of `cx` is absent, or when
    a mitochondrial gene's reference chain is not that gene."""
    d = MMCIF2Dict(str(cif))
    unp = {}
    for ent, db, code in zip(d["_struct_ref.entity_id"], d["_struct_ref.db_name"], d["_struct_ref.db_code"]):
        if db == "UNP":
            unp.setdefault(ent, set()).add(code)
    out = {}
    for ent, chains in zip(d["_entity_poly.entity_id"], d["_entity_poly.pdbx_strand_id"]):
        codes = unp.get(ent, set())
        if len(codes) != 1 or next(iter(codes)) not in cx.uniprot:
            raise SystemExit(f"{cif.name}: entity {ent} (chains {chains}) has UniProt reference {sorted(codes) or 'none'}, "
                             f"which complex {cx.name}'s table does not name; refusing to guess its gene")
        for c in chains.split(","):
            out[c.strip()] = cx.uniprot[next(iter(codes))]
    absent = sorted(set(cx.uniprot.values()) - set(out.values()))
    if absent:
        raise SystemExit(f"{cif.name}: no chain of {', '.join(absent)}, which complex {cx.name}'s table expects")
    for g, c in cx.mt.items():
        if out.get(c) != g:
            raise SystemExit(f"{cif.name}: chain {c} is {out.get(c)}, not {g} as complex {cx.name}'s table says")
    if set(cx.nuclear) != {g for g in out.values() if g not in MT_GENES} or set(cx.mt) != {g for g in out.values() if g in MT_GENES}:
        raise SystemExit(f"{cif.name}: its genes are not complex {cx.name}'s nuclear and mitochondrial lists")
    return out
