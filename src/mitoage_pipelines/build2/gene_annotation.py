"""Protein genes a record leaves out, and what build 2 does about each: the hand-checked lists of the site's
pipeline/sources/gene-annotation.ts (its decision D17, commit 03c0949). Like POSSIBLE_DUPLICATES, these are
curated inputs of the build, not rules, so they are copied here as data; the code that applies them is written
from the rules (genome.apply_gene_patches, assemble.gene_flags).

- ANNOTATION_PATCHES: a gene intact in a record's sequence but not annotated. The build adds it to that exact
  accession.version only, and refuses the patch if the record annotates the gene itself or if the location does
  not read as one open reading frame.
- GENE_GAPS, LINEAGE_GENE_GAPS: a gene the build leaves out, with the reason, published as a review flag. A live
  species missing a gene that neither list names is flagged gene_unexplained_absence and listed for review.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class AnnotationPatch:
    accession: str
    gene: str
    location: str      # GenBank location syntax, in the record's coordinates
    evidence: str      # how the location was found and checked; the changelog row's reason


@dataclass(frozen=True)
class GeneGap:
    accession: str
    gene: str
    reason: str        # absent_from_record, frameshifted_in_record, not_annotated or lost_in_lineage
    detail: str


@dataclass(frozen=True)
class LineageGap:
    gene: str
    rank: str          # "phylum" or "class"
    taxon: str
    reason: str
    detail: str


ANNOTATION_PATCHES = (
    AnnotationPatch(
        "NC_027847.1", "ND6", "complement(17376..17891)",
        "Psittacus erithacus: the record annotates no ND6. An open reading frame of 171 codons and a TAA stop on the "
        "minus strand, between tRNA-Pro (complement(17303..17373)) and tRNA-Glu (complement(17893..17962)), the avian "
        "position of ND6, translates to a protein matching Amazona aestiva ND6 over its whole length; found by a "
        "six-frame search with the Amazona protein, 2026-11-15"),
)

GENE_GAPS = (
    GeneGap("NC_004815.1", "ND5", "absent_from_record",
            "Sphenodon punctatus: the record's comment says the genome lacks ND5, tRNA-His and tRNA-Thr (Rest et al. 2003)"),
    GeneGap("NC_010689.1", "ND6", "absent_from_record",
            "Chionodraco myersi: the record says the genome lacks ND6 and tRNA-Glu; in Antarctic notothenioids ND6 moved into "
            "the control region (Zhuang & Cheng 2010, doi:10.1093/molbev/msq026), which this record does not carry it in"),
    GeneGap("NC_020020.1", "ND2", "absent_from_record",
            "Homarus gammarus: a provisional record of 14,316 bp that starts inside COX1; ND2 is not in its sequence"),
    GeneGap("NC_033966.1", "CYTB", "frameshifted_in_record",
            "Pluvialis fulva: misc_feature 13637..14774 'nonfunctional CytB due to frameshift'"),
    GeneGap("NC_008140.1", "ND1", "frameshifted_in_record",
            "Podiceps cristatus: misc_feature 2716..3694 'similar to NADH dehydrogenase subunit 1', reading in two frames"),
)

LINEAGE_GENE_GAPS = (
    LineageGap("ATP8", "phylum", "Nematoda", "lost_in_lineage",
               "nematode mtDNA has no ATP8 (Okimoto et al. 1992, doi:10.1093/genetics/130.3.471)"),
    LineageGap("ATP8", "class", "Bivalvia", "not_annotated",
               "many bivalve records annotate no ATP8; whether the gene is lost or an unannotated short reading frame was not checked"),
)
