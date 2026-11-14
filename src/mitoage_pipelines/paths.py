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
EMBEDDINGS = INPUT / "embeddings"         # ESM-2 mean embeddings of both partners (prepared data only)
NCBI = INPUT / "ncbi"                     # NCBI ortholog sets, one folder per gene: ncbi/<GENE>/ncbi_dataset/data
MANIFEST = INPUT / "MANIFEST.json"        # what was downloaded, from where and when

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
B2_EURING = B2_IN / "euring"              # status.json: what the attempt to download it found
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
MITOPORT = INTERIM / "mitoport"
MITOPORT_SEL = MITOPORT / "sel"           # one protein (.faa) and coding sequence (.fna) per species and gene
MITOPORT_MAFFT = MITOPORT / "mafft"       # MAFFT's alignments as it wrote them
MITOPORT_ALN = MITOPORT / "aln"           # filtered alignments on the complete species matrix
MITOPORT_BL = MITOPORT / "bl"             # IQ-TREE branch lengths per protein on the fixed topology
MITOPORT_TOPOLOGY = MITOPORT / "topology.nwk"
MITOPORT_GENES = MITOPORT / "genes.txt"
MITOPORT_SPECIES = MITOPORT / "species.txt"

B2_INTERIM = INTERIM / "build2"           # refseq_records.parquet, observations.parquet, matches.parquet, plan.json, ...

# --- output: result tables --------------------------------------------------------------------
LONGEVITY_OUT = OUTPUT / "longevity"      # composition.tsv controls.tsv fourfold.tsv controls_composition.tsv
                                          # sites.tsv pgls.<Class>.tsv aa_pgls.<Class>.tsv agreement.tsv columns.<GENE>.tsv
RELAX_OUT = LONGEVITY_OUT / "relax"       # <GENE>.relax.json, as HyPhy writes them
PEPTIDES_OUT = OUTPUT / "peptides"        # <peptide>.tsv, summary.txt
MITOPORT_OUT = OUTPUT / "mitoport"        # erc.tsv


def tree_file(cls: str) -> Path:
    """The IQ-TREE tree of one class, as `mitoage longevity tree` writes it."""
    return TREE / f"{cls}.iq.treefile"


def outgroup_file(cls: str) -> Path:
    return TREE / f"{cls}.outgroup.txt"


def aa_pgls_table(cls: str) -> Path:
    """H1 on the dated tree, one class, as `mitoage longevity aa-pgls` writes it."""
    return LONGEVITY_OUT / f"aa_pgls.{cls}.tsv"


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

