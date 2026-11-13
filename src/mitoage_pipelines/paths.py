"""Where every file lives. The only module that knows the folder layout.

    data/input     downloads, never edited: the MitoAge build, its dated species tree, the NCBI ortholog sets
    data/interim   made here on the way to a result: alignments, trees, selections, inputs of other programs
    data/output    result tables

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
NCBI = INPUT / "ncbi"                     # NCBI ortholog sets, one folder per gene: ncbi/<GENE>/ncbi_dataset/data
MANIFEST = INPUT / "MANIFEST.json"        # what was downloaded, from where and when

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

# --- output: result tables --------------------------------------------------------------------
LONGEVITY_OUT = OUTPUT / "longevity"      # composition.tsv controls.tsv fourfold.tsv controls_composition.tsv
                                          # sites.tsv pgls.<Class>.tsv agreement.tsv columns.<GENE>.tsv
RELAX_OUT = LONGEVITY_OUT / "relax"       # <GENE>.relax.json, as HyPhy writes them
PEPTIDES_OUT = OUTPUT / "peptides"        # <peptide>.tsv, summary.txt
MITOPORT_OUT = OUTPUT / "mitoport"        # erc.tsv


def tree_file(cls: str) -> Path:
    """The IQ-TREE tree of one class, as `mitoage longevity tree` writes it."""
    return TREE / f"{cls}.iq.treefile"


def outgroup_file(cls: str) -> Path:
    return TREE / f"{cls}.outgroup.txt"


def pgls_table(cls: str, gene_tree: bool = False) -> Path:
    return LONGEVITY_OUT / (f"pgls.gene-tree.{cls}.tsv" if gene_tree else f"pgls.{cls}.tsv")
