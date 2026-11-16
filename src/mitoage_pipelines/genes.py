"""Gene lists, in one place."""

# the 13 proteins of the vertebrate mitochondrial genome, named as the build's alignment files are
MT_GENES = "ATP6 ATP8 COX1 COX2 COX3 CYTB ND1 ND2 ND3 ND4 ND4L ND5 ND6".split()
# ND6 is on the other strand: base composition at fourfold sites is read from the other twelve
MT_HEAVY_STRAND = [g for g in MT_GENES if g != "ND6"]

COMPLEX_IV_NUCLEAR = "COX4I1 COX5A COX5B COX6A1 COX6B1 COX6C COX7A2 COX7B COX7C COX8A NDUFA4".split()
COMPLEX_III_NUCLEAR = "UQCRC1 UQCRC2 CYC1 UQCRFS1 UQCRB UQCRQ UQCRH UQCR10 UQCR11".split()
CYTOCHROME_C = ["CYCS"]
COMPLEX_II = "SDHA SDHB SDHC SDHD".split()
MATRIX = "CS MDH2 SOD2 HSPD1 ACO2 IDH2 FH".split()
CYTOSOLIC_CONTROLS = ("CSTB FABP3 MB LYZ RBP4 TTR B2M CD59 PFN1 CFL1 GSTP1 PRDX1 PRDX2 TXN PPIA FKBP1A NME1 PEBP1 "
                      "LGALS1 CRYAB HSPB1 TPT1 PGAM1 CA2 GPX1 SOD1 TPI1 LDHA PGK1 GAPDH").split()

# the groups of the rate-covariation test, by the names its table carries. The cytosolic controls are
# not listed here: the test takes as control every nuclear protein in its matrix that is in no group.
GROUPS = {
    "complex IV nuclear": COMPLEX_IV_NUCLEAR,
    "complex III nuclear": COMPLEX_III_NUCLEAR,
    "cytochrome c": CYTOCHROME_C,
    "complex II (no mt subunit)": COMPLEX_II,
    "matrix, not respiratory": MATRIX,
}

# every nuclear gene whose ortholog set `mitoage download orthologs` fetches
NUCLEAR = COMPLEX_IV_NUCLEAR + COMPLEX_III_NUCLEAR + CYTOCHROME_C + COMPLEX_II + MATRIX + CYTOSOLIC_CONTROLS

# the nuclear subunits of the solved human complexes I (5XTD) and V (8H9S), as `mitoport complexes` names their
# entities; not among the genes H5 downloads by default, fetched by name for the contact maps of complexes I and V
COMPLEX_I_NUCLEAR = ("NDUFV1 NDUFV2 NDUFV3 NDUFS1 NDUFS2 NDUFS3 NDUFS4 NDUFS5 NDUFS6 NDUFS7 NDUFS8 "
                     "NDUFA1 NDUFA2 NDUFA3 NDUFA5 NDUFA6 NDUFA7 NDUFA8 NDUFA9 NDUFA10 NDUFA11 NDUFA12 NDUFA13 NDUFAB1 "
                     "NDUFB1 NDUFB2 NDUFB3 NDUFB4 NDUFB5 NDUFB6 NDUFB7 NDUFB8 NDUFB9 NDUFB10 NDUFB11 NDUFC1 NDUFC2").split()
COMPLEX_V_NUCLEAR = ("ATP5F1A ATP5F1B ATP5F1C ATP5F1D ATP5F1E ATP5PB ATP5PD ATP5PF ATP5PO ATP5MC1 ATP5ME ATP5MF "
                     "ATP5MG ATP5MJ ATP5IF1").split()

# H19's membrane control: integral proteins of the inner mitochondrial membrane that are in no respiratory complex and
# have no mtDNA-encoded partner (carriers, the pyruvate carrier, protein import, MICOS, LETM1, the calcium uniporter,
# prohibitins). TIMM23 is left out: NCBI's symbol lookup returns a gene with 134 vertebrate orthologs.
INNER_MEMBRANE = ("SLC25A1 SLC25A11 SLC25A12 SLC25A20 SLC25A22 MPC1 MPC2 TIMM17A TIMM22 TIMM50 IMMT LETM1 MCU "
                  "PHB1 PHB2").split()

# every nuclear gene H19 reads from NCBI's vertebrate ortholog sets (`mitoport classes`)
H19_NUCLEAR = (COMPLEX_I_NUCLEAR + COMPLEX_III_NUCLEAR + COMPLEX_IV_NUCLEAR + COMPLEX_V_NUCLEAR + CYTOCHROME_C
               + COMPLEX_II + MATRIX + INNER_MEMBRANE + CYTOSOLIC_CONTROLS)

# mitochondrial subunits of the complexes with nuclear partners above, named as in the rate matrix
MT_OF_COMPLEX = {"IV": ["MT-COX1", "MT-COX2", "MT-COX3"], "III": ["MT-CYTB"]}
