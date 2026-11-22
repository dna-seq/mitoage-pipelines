"""Where every file lives. The only module that knows the folder layout.

    data/input     downloads, never edited: the MitoAge build, its dated species tree, the NCBI ortholog sets
    data/interim   made here on the way to a result: alignments, trees, selections, inputs of other programs
    data/output    result tables
    data/reference our own result tables, from `download prepared`, to compare a rerun with

The root is `data/` in the repository (or in the working directory, when the package is installed
outside a checkout); the environment variable MITOAGE_PIPELINES_DATA overrides it.
"""
import os
from pathlib import Path

ENV = "MITOAGE_PIPELINES_DATA"


def data_root() -> Path:
    env = os.environ.get(ENV)
    if env:
        return Path(env).expanduser().resolve()
    repo = Path(__file__).resolve().parents[2]          # src/mitoage_pipelines/paths.py -> the checkout
    if (repo / "pyproject.toml").is_file():
        return repo / "data"
    return Path.cwd() / "data"


ROOT = data_root()

INPUT = ROOT / "input"
INTERIM = ROOT / "interim"
OUTPUT = ROOT / "output"

# --- input: downloads -------------------------------------------------------------------------
BUILD = INPUT / "mitoage"                 # the Hugging Face dataset longevity-genie/mitoage
SITE_TREE = INPUT / "mitoage-site" / "phylogeny.json"   # the build's dated species tree, from the site's repository
EMBEDDINGS = INPUT / "embeddings"         # ESM-2 embeddings of both partners, mean and per residue (prepared data only)
NCBI = INPUT / "ncbi"                     # NCBI ortholog sets, one folder per gene: ncbi/<GENE>/ncbi_dataset/data
NCBI_VERT = INPUT / "ncbi-vertebrata"     # the same for complex IV's nuclear subunits, within vertebrates (H12)
MANIFEST = INPUT / "MANIFEST.json"        # what was downloaded, from where and when
PDB = INPUT / "pdb"                       # solved structures from RCSB, <ID>.cif as served
CYBRID_GB = INPUT / "cybrid"              # records.gb: the mitogenomes of the cybrid hosts and donors (H16)
STABDDG = INPUT / "stabddg"              # the three weight files of StaB-ddG at its pinned commit (H28)
TECH_GB = INPUT / "technology"            # refseq/ and source/: GenBank records the sequencing technology is read from (H8's control)
HIBERNATION = INPUT / "hibernation"        # the bat table of Wilkinson & Adams 2019 (DRUM) and Heldstab et al. 2018 (Dryad, placed by hand)

# --- input: the sources of build 2 (docs/BUILD2_REPRODUCTION.md) --------------------------------
B2_IN = INPUT / "build2"
B2_HF = B2_IN / "hf"                      # longevity-genie/mitoage at the pinned revision: build1/ and the build 2 tables
B2_BUILD1 = B2_HF / "build1"              # MitoAge build 1, the 2015 CSV export as parquet
B2_REFSEQ = B2_IN / "refseq"              # release<N>/mitochondrion.1.genomic.gbff.gz, or efetch/records.gb.gz
B2_TAXDUMP = B2_IN / "taxdump"            # <date>/nodes.dmp names.dmp merged.dmp rankedlineage.dmp
B2_ANAGE = B2_IN / "anage" / "dataset.zip"
B2_FISHBASE = B2_IN / "fishbase"          # <version>/<table>.parquet
B2_SEALIFEBASE = B2_IN / "sealifebase"    # <version>/<table>.parquet
B2_USGS = B2_IN / "usgs_bbl"              # Longevity.php.html and longevity.csv, the table read from it
B2_EURING = B2_IN / "euring"              # EURING_longevity_list_20230901.pdf (the Internet Archive's copy)
B2_GENBANK = B2_IN / "genbank"            # cited/ (build 1's records), contig/ (sequence of CONTIG-only records), primary/ (INSDC)
B2_LEGACY = B2_IN / "legacy"              # pages of the old site with its species ids
B2_TREE_IN = B2_IN / "tree"               # Open Tree of Life and TimeTree answers, one JSON file per request

# --- interim ----------------------------------------------------------------------------------
ALIGNERS = INTERIM / "aligners"           # <GENE>.<tool>.afa: one gene aligned by several tools
TREE = INTERIM / "tree"                   # dated.nwk, quality_ok.txt; <Class>.faa, <Class>.outgroup.txt, <Class>.iq.treefile
DATED_TREE = TREE / "dated.nwk"           # the dated species tree as newick, tips named by species slug
QUALITY_OK = TREE / "quality_ok.txt"      # species whose lifespan record is of high or acceptable quality
RELAX = INTERIM / "relax"                 # species.tsv, tree.nwk, <GENE>.fna: the inputs of HyPhy RELAX
COEVOL = INTERIM / "coevol"               # mammals.ali, mammals.tree, mammals.lht
COEVOL_FIXED = INTERIM / "coevol-fixed"   # the same, with the dated tree's ages kept as branch lengths (-fixtimes)
MITOPORT = INTERIM / "mitoport"
MITOPORT_SEL = MITOPORT / "sel"           # one protein (.faa) and coding sequence (.fna) per species and gene
MITOPORT_MAFFT = MITOPORT / "mafft"       # MAFFT's alignments as it wrote them
MITOPORT_ALN = MITOPORT / "aln"           # filtered alignments on the complete species matrix
MITOPORT_BL = MITOPORT / "bl"             # IQ-TREE branch lengths per protein on the fixed topology
MITOPORT_TOPOLOGY = MITOPORT / "topology.nwk"
MITOPORT_GENES = MITOPORT / "genes.txt"
MITOPORT_SPECIES = MITOPORT / "species.txt"
MITOPORT_IFACE_ALN = MITOPORT / "iface_aln"   # every nuclear selection, MAFFT L-INS-i --anysymbol, unfiltered
PARTNERS = MITOPORT / "partners"          # sel/ and aln/: complex IV's nuclear subunits in four classes (H12)
MITOPORT_IFACE = MITOPORT / "iface"       # contacts.tsv, sites.tsv, control_pairs.tsv, chains.tsv, species.txt: the complex IV contact map
STRUCTURE = INTERIM / "structure"
STRUCTURE_PDBS = STRUCTURE / "pdbs"       # <GENE>/<species>.pdb predicted by ESMFold2, plddt.tsv
STRUCTURE_CHAINS = STRUCTURE / "chains"   # <GENE>.pdb: one chain of a solved human structure per gene
STRUCTURE_TMP = STRUCTURE / "msa"         # the sequence and structural alignments the check compares
EMBED_COX = STRUCTURE / "embed_cox.parquet"   # ESM-2 means over complex IV's contact and control residues (H12)

AUDIT = INTERIM / "audit"                 # records.parquet: the RefSeq record of every species of the build (H13)
B2_INTERIM = INTERIM / "build2"           # refseq_records.parquet, observations.parquet, matches.parquet, plan.json, ...

# --- output: result tables --------------------------------------------------------------------
LONGEVITY_OUT = OUTPUT / "longevity"      # composition.tsv controls.tsv fourfold.tsv controls_composition.tsv
                                          # sites.tsv pgls.<Class>.tsv pgls_power.<Class>.tsv aa_pgls.<Class>.tsv
                                          # agreement.tsv columns.<GENE>.tsv
RELAX_OUT = LONGEVITY_OUT / "relax"       # <GENE>.relax.json, as HyPhy writes them
PEPTIDES_OUT = OUTPUT / "peptides"        # <peptide>.tsv, summary.txt
MITOPORT_OUT = OUTPUT / "mitoport"        # erc.tsv, coembed.tsv
IFACE_OUT = MITOPORT_OUT / "iface"        # coupling.tsv, lifespan.tsv
CYBRID_OUT = MITOPORT_OUT / "cybrid"      # positions.tsv, scores.tsv, calibration.tsv, summary.tsv (H16); coev_*, esm_*, h17_* (H17)
STRUCTURE_OUT = OUTPUT / "structure"      # within_species.tsv, chains.tsv, structural_msa.tsv
ADAPT_OUT = OUTPUT / "adaptations"        # traits, pairs, fits, selection/, coupling/, readout (H21-H25)
ADAPT_LONGEVITY_OUT = OUTPUT / "adapt-longevity"   # the same steps for H26's trait, long-lived for body mass
AUDIT_OUT = OUTPUT / "audit"              # inputs, strand, order, code, pairs, nearest (.tsv): the audit of the records (H13)
COEVOL_OUT = LONGEVITY_OUT / "coevol"     # chain<N>.* as Coevol writes them, chain<N>.cov from readcoevol
COEVOL_FIXED_OUT = LONGEVITY_OUT / "coevol-fixed"   # the chains run with the ages held (-fixtimes)
THERMAL_OUT = LONGEVITY_OUT / "thermal"      # H20: traits, pairs, pair_fits, composition, lehmann, readout
BRANCH_RATES_OUT = LONGEVITY_OUT / "branch_rates"   # H3's two-step check: fit.json, branches.tsv, pgls.tsv


def tree_file(cls: str) -> Path:
    """The IQ-TREE tree of one class, as `mitoage longevity tree` writes it."""
    return TREE / f"{cls}.iq.treefile"


def outgroup_file(cls: str) -> Path:
    return TREE / f"{cls}.outgroup.txt"


def aa_pgls_table(cls: str) -> Path:
    """H1 on the dated tree, one class, as `mitoage longevity aa-pgls` writes it."""
    return LONGEVITY_OUT / f"aa_pgls.{cls}.tsv"


def pgls_power_table(cls: str) -> Path:
    """The power of the column test of one class, as `mitoage longevity pgls-power` writes it."""
    return LONGEVITY_OUT / f"pgls_power.{cls}.tsv"


def pgls_table(cls: str, gene_tree: bool = False) -> Path:
    return LONGEVITY_OUT / (f"pgls.gene-tree.{cls}.tsv" if gene_tree else f"pgls.{cls}.tsv")


# --- output: the reproduction of build 2 --------------------------------------------------------
B2_OUT = OUTPUT / "build2"                # species, region_stats, codon_usage, genome_features, lifespan_observations,
                                          # changelog (.parquet), review_list.csv, phylogeny.json, alignments/, compare/


def refseq_release(release: str) -> Path:
    """The RefSeq mitochondrion release file of one release, as `download refseq` writes it."""
    return B2_REFSEQ / f"release{release}" / "mitochondrion.1.genomic.gbff.gz"


def taxdump(label: str) -> Path:
    """The folder of one NCBI taxonomy dump (an archive date such as 2026-11-01, or current-<date>)."""
    return B2_TAXDUMP / label

