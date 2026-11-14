"""The `mitoage` command: every step, with the folders of `paths.py` as defaults.

Order of a full run:
  mitoage download build && mitoage download orthologs
  mitoage longevity controls && mitoage longevity fourfold && mitoage longevity composition
  mitoage longevity sites
  mitoage download tree && mitoage longevity dated-tree && mitoage longevity pgls
  mitoage longevity relax-prep && mitoage longevity relax
  mitoage longevity coevol-prep
  mitoage peptides pilot
  mitoage mitoport select && mitoage mitoport align && mitoage mitoport branch-lengths && mitoage mitoport erc
  mitoage download prepared --topic embeddings && mitoage mitoport coembed
Heavy imports happen inside the commands, so that `--help` is quick.
"""
from pathlib import Path
from typing import Annotated, Optional

import typer

from . import genes as gene_lists
from . import paths

app = typer.Typer(no_args_is_help=True, pretty_exceptions_enable=False,
                  help="MitoAge analyses: download the data, then run each step.")
download_app = typer.Typer(no_args_is_help=True, help="Download the inputs into data/input.")
longevity_app = typer.Typer(no_args_is_help=True, help="Longevity in the alignments of the 13 mitochondrial proteins.")
peptides_app = typer.Typer(no_args_is_help=True, help="Peptides encoded inside mtDNA.")
mitoport_app = typer.Typer(no_args_is_help=True, help="Rate covariation of mitochondrial proteins with their nuclear partners.")
app.add_typer(download_app, name="download")
app.add_typer(longevity_app, name="longevity")
app.add_typer(peptides_app, name="peptides")
app.add_typer(mitoport_app, name="mitoport")

Build = Annotated[Path, typer.Option("--build", help="the MitoAge build (species.parquet, alignments/)")]
Out = Annotated[Path, typer.Option("--out", help="folder the tables are written to")]
Composition = Annotated[Path, typer.Option("--composition", help="composition.tsv, written by `longevity controls`")]
Fourfold = Annotated[Path, typer.Option("--fourfold", help="fourfold.tsv, written by `longevity fourfold`")]
Class = Annotated[str, typer.Option("--class", help="taxonomic class, as in species.parquet")]
AltDir = Annotated[Optional[Path], typer.Option("--alignments-dir", help="read <dir>/<GENE>.<tool>.afa instead of the build's protein alignments (with --tool)")]
AltTool = Annotated[Optional[str], typer.Option("--tool", help="the tool's name in those file names")]

COMPOSITION = paths.LONGEVITY_OUT / "composition.tsv"
FOURFOLD = paths.LONGEVITY_OUT / "fourfold.tsv"


# --- download ---------------------------------------------------------------------------------
@download_app.command("build")
def download_build(
    dest: Annotated[Path, typer.Option(help="where the dataset goes")] = paths.BUILD,
    manifest: Annotated[Path, typer.Option(help="the download record")] = paths.MANIFEST,
    revision: Annotated[Optional[str], typer.Option(help="a commit of the dataset; the newest when omitted")] = None,
) -> None:
    """The MitoAge build from Hugging Face (longevity-genie/mitoage)."""
    from . import download
    download.build(dest, manifest, revision)


@download_app.command("prepared")
def download_prepared(
    topic: Annotated[str, typer.Option(help="partners, embeddings, aligners, results, or all")] = "all",
    revision: Annotated[Optional[str], typer.Option(help="a commit of the dataset, to pin what is read (one topic only)")] = None,
) -> None:
    """What we fetched and computed, from Hugging Face, by topic: start any step without redoing those before it.
    Our result tables go to data/reference, for comparison with a rerun's data/output."""
    from . import download
    topics = list(download.PREPARED) if topic == "all" else [topic]
    if topic != "all" and topic not in download.PREPARED:
        raise typer.BadParameter(f"unknown topic {topic!r}; one of {', '.join(download.PREPARED)} or all")
    for t in topics:
        download.prepared(t, paths.ROOT, paths.MANIFEST, revision if len(topics) == 1 else None)


@download_app.command("tree")
def download_tree(ref: Annotated[str, typer.Option(help="branch, tag or commit of the site's repository")] = "main") -> None:
    """The build's dated species tree (phylogeny.json) from the site's repository into data/input/mitoage-site."""
    from . import download
    download.site_tree(paths.SITE_TREE, paths.MANIFEST, ref)


@download_app.command("orthologs")
def download_orthologs(
    gene: Annotated[Optional[list[str]], typer.Option("--gene", help="a gene symbol; repeat for several. All the nuclear genes of the analysis when omitted")] = None,
    dest: Annotated[Path, typer.Option(help="where the ortholog sets go, one folder per gene")] = paths.NCBI,
    manifest: Annotated[Path, typer.Option(help="the download record")] = paths.MANIFEST,
    force: Annotated[bool, typer.Option("--force", help="download a gene again even when it is there")] = False,
) -> None:
    """Mammalian ortholog sets (protein and CDS) of the nuclear genes, from NCBI Datasets."""
    from . import download
    download.orthologs(gene or gene_lists.NUCLEAR, dest, manifest, skip_existing=not force)


# --- longevity --------------------------------------------------------------------------------
@longevity_app.command("compare-aligners")
def longevity_compare_aligners(
    files: Annotated[list[Path], typer.Argument(help="alignments of the same sequences, named <GENE>.<tool>.afa")],
    out: Out = paths.LONGEVITY_OUT,
    ref: Annotated[str, typer.Option(help="the tool the others are compared with, column by column")] = "mafft_fast",
) -> None:
    """Agreement between alignments of one gene made by different aligners."""
    from .longevity import compare
    compare.main(files, out, ref)


@longevity_app.command("controls")
def longevity_controls(data: Build = paths.BUILD, out: Out = paths.LONGEVITY_OUT) -> None:
    """Amino-acid content of the 13 proteins against maximum lifespan: plain, within families, mixed model."""
    from .longevity import controls
    controls.main(data, out)


@longevity_app.command("fourfold")
def longevity_fourfold(data: Build = paths.BUILD, out: Out = paths.LONGEVITY_OUT) -> None:
    """Base composition at fourfold-degenerate sites per species."""
    from .longevity import fourfold
    fourfold.main(data, out)


@longevity_app.command("composition")
def longevity_composition(data: Build = paths.BUILD, composition: Composition = COMPOSITION, out: Out = paths.LONGEVITY_OUT) -> None:
    """The amino-acid tests in mammals with fourfold base composition in the model."""
    from .longevity import composition as step
    step.main(data, composition, out)


@longevity_app.command("sites")
def longevity_sites(
    data: Build = paths.BUILD, composition: Composition = COMPOSITION, fourfold: Fourfold = FOURFOLD,
    out: Annotated[Path, typer.Option("--out", help="the table to write")] = paths.LONGEVITY_OUT / "sites.tsv",
    alignments_dir: AltDir = None, tool: AltTool = None,
) -> None:
    """Every alignment column against lifespan within families, no tree (the cautionary version: inflated)."""
    from .longevity import sites
    sites.main(data, composition, fourfold, out, alignments_dir, tool)


@longevity_app.command("concat")
def longevity_concat(
    data: Build = paths.BUILD, composition: Composition = COMPOSITION,
    out: Annotated[Path, typer.Option("--out", help="folder for <Class>.faa and <Class>.outgroup.txt")] = paths.TREE,
    cls: Class = "Mammalia",
) -> None:
    """Concatenated proteins of one class plus an outgroup, for the tree."""
    from .longevity import concat
    concat.main(data, composition, out, cls)


@longevity_app.command("tree")
def longevity_tree(
    tree_dir: Annotated[Path, typer.Option("--tree-dir", help="folder holding <Class>.faa; the tree is written next to it")] = paths.TREE,
    cls: Class = "Mammalia",
    model: Annotated[Optional[str], typer.Option(help="IQ-TREE model; mtMAM+F+G4 for Mammalia and mtVer+F+G4 for Aves when omitted")] = None,
    threads: Annotated[str, typer.Option(help="IQ-TREE -T")] = "AUTO",
    seed: Annotated[Optional[int], typer.Option(help="IQ-TREE --seed; a random one when omitted, as in the published run")] = None,
) -> None:
    """The tree of one class with IQ-TREE (from PATH): -m <model> -fast."""
    from .longevity import tree
    print(tree.main(tree_dir, cls, model, threads, seed))


@longevity_app.command("dated-tree")
def longevity_dated_tree(data: Build = paths.BUILD,
                         site_tree: Annotated[Path, typer.Option(help="phylogeny.json, from `download tree`")] = paths.SITE_TREE) -> None:
    """The dated species tree as newick, and the species with a lifespan record of high or acceptable quality."""
    from .longevity import dated_tree
    dated_tree.main(data, site_tree, paths.DATED_TREE, paths.QUALITY_OK)


@longevity_app.command("pgls")
def longevity_pgls(
    data: Build = paths.BUILD, composition: Composition = COMPOSITION, fourfold: Fourfold = FOURFOLD,
    cls: Class = "Mammalia",
    gene_tree: Annotated[bool, typer.Option("--gene-tree", help="use the tree built from the 13 proteins (`longevity tree`) and every lifespan record: NOT the analysis as fixed, kept to show the result's dependence on the tree")] = False,
    tree: Annotated[Optional[Path], typer.Option(help="the tree; the dated species tree when omitted")] = None,
    outgroup: Annotated[Optional[Path], typer.Option(help="with --gene-tree: the outgroup's species")] = None,
    out: Annotated[Optional[Path], typer.Option("--out", help="the table to write; <output>/longevity/pgls.<Class>.tsv when omitted")] = None,
    alignments_dir: AltDir = None, tool: AltTool = None,
) -> None:
    """The column test as a phylogenetic regression (Pagel's lambda) on the dated species tree."""
    from .longevity import pgls
    if gene_tree:
        pgls.main(data, composition, fourfold, tree or paths.tree_file(cls), outgroup or paths.outgroup_file(cls),
                  out or paths.pgls_table(cls, gene_tree=True), cls, alignments_dir, tool)
    else:
        pgls.main(data, composition, fourfold, tree or paths.DATED_TREE, None, out or paths.pgls_table(cls),
                  cls, alignments_dir, tool, quality_file=paths.QUALITY_OK)


@longevity_app.command("relax-prep")
def longevity_relax_prep(
    data: Build = paths.BUILD, composition: Composition = COMPOSITION,
    tree: Annotated[Path, typer.Option(help="the species tree; the dated one, from `longevity dated-tree`")] = paths.DATED_TREE,
    out: Annotated[Path, typer.Option("--out", help="folder for species.tsv, tree.nwk and <GENE>.fna")] = paths.RELAX,
    n: Annotated[int, typer.Option(help="species per order at each end, at most")] = 10,
    per_family: Annotated[int, typer.Option(help="species of one family at each end, at most")] = 3,
    exclude: Annotated[Optional[Path], typer.Option(help="species.tsv of an earlier draw: its species are left out")] = None,
) -> None:
    """Long-lived against short-lived mammals, balanced within orders: tree and codon files for HyPhy RELAX.

    The second, disjoint draw of the analysis is the same command with --exclude <first draw>/species.tsv
    and another --out."""
    from .longevity import relax_prep
    dated = tree == paths.DATED_TREE
    relax_prep.main(data, composition, tree, out, n, per_family, exclude,
                    quality_file=paths.QUALITY_OK if dated else None, in_tree_only=dated)


@longevity_app.command("relax")
def longevity_relax(
    relax_dir: Annotated[Path, typer.Option("--relax-dir", help="the folder `relax-prep` wrote")] = paths.RELAX,
    out: Annotated[Path, typer.Option("--out", help="folder for <GENE>.relax.json")] = paths.RELAX_OUT,
    gene: Annotated[Optional[list[str]], typer.Option("--gene", help="a gene; repeat for several. All 13 when omitted")] = None,
    cpus: Annotated[int, typer.Option(help="HyPhy CPU=")] = 1,
    minimal: Annotated[bool, typer.Option("--minimal", help="--models Minimal: the test alone")] = False,
) -> None:
    """HyPhy RELAX per gene (hyphy from PATH): --code Vertebrate-mtDNA --test Foreground --reference Reference."""
    from .longevity import relax
    relax.main(relax_dir, out, gene, cpus, minimal)


@longevity_app.command("coevol-prep")
def longevity_coevol_prep(
    data: Build = paths.BUILD, composition: Composition = COMPOSITION,
    tree: Annotated[Path, typer.Option(help="the species tree; the dated one, from `longevity dated-tree`")] = paths.DATED_TREE,
    outgroup: Annotated[Optional[Path], typer.Option(help="only with a tree built from the genes: its outgroup's species, to root it")] = None,
    out: Annotated[Path, typer.Option("--out", help="folder for mammals.ali, mammals.tree, mammals.lht")] = paths.COEVOL,
    per_family: Annotated[int, typer.Option(help="species of one family, at most")] = 1,
) -> None:
    """Inputs for Coevol: codons, tree, lifespan and mass."""
    from .longevity import coevol_prep
    coevol_prep.main(data, composition, tree, outgroup, out, per_family,
                     quality_file=paths.QUALITY_OK if outgroup is None else None)


@longevity_app.command("carriers")
def longevity_carriers(data: Build = paths.BUILD, composition: Composition = COMPOSITION) -> None:
    """Which families carry the minor residue at ND5 column 316 and CYTB column 83."""
    from .longevity import carriers
    carriers.main(data, composition)


# --- peptides ---------------------------------------------------------------------------------
@peptides_app.command("pilot")
def peptides_pilot(data: Build = paths.BUILD, out: Out = paths.PEPTIDES_OUT) -> None:
    """Are the reading frames of published mtDNA-encoded peptides kept across species?"""
    from .peptides import pilot
    pilot.main(data, out)


# --- mitoport ---------------------------------------------------------------------------------
Ncbi = Annotated[Path, typer.Option("--ncbi", help="the NCBI ortholog sets, one folder per gene")]
Sel = Annotated[Path, typer.Option("--sel", help="folder of the selections: <GENE>.faa, <GENE>.fna")]
AlnDir = Annotated[Path, typer.Option("--aln-dir", help="folder of the filtered alignments")]
Topology = Annotated[Path, typer.Option(help="the fixed topology")]
GenesFile = Annotated[Path, typer.Option("--genes", help="the proteins of the matrix, one per line")]
BlDir = Annotated[Path, typer.Option("--bl-dir", help="folder of IQ-TREE's branch lengths per protein")]


@mitoport_app.command("overlap")
def mitoport_overlap(data: Build = paths.BUILD, ncbi: Ncbi = paths.NCBI) -> None:
    """How many of the build's mammals each NCBI ortholog set covers."""
    from .mitoport import overlap
    overlap.main(data, ncbi)


@mitoport_app.command("select")
def mitoport_select(
    data: Build = paths.BUILD, ncbi: Ncbi = paths.NCBI,
    out: Annotated[Path, typer.Option("--out", help="folder for <GENE>.faa and <GENE>.fna")] = paths.MITOPORT_SEL,
) -> None:
    """One protein and coding sequence per species and gene from the NCBI ortholog sets."""
    from .mitoport import select
    select.main(data, ncbi, out)


@mitoport_app.command("align")
def mitoport_align(
    data: Build = paths.BUILD, sel: Sel = paths.MITOPORT_SEL,
    tree: Annotated[Path, typer.Option(help="the species tree; the dated one, from `longevity dated-tree`")] = paths.DATED_TREE,
    outgroup: Annotated[Optional[Path], typer.Option(help="only with a tree built from the genes: its outgroup's species, to root it")] = None,
    mafft_dir: Annotated[Path, typer.Option("--mafft-dir", help="folder for MAFFT's alignments as written")] = paths.MITOPORT_MAFFT,
    aln_dir: AlnDir = paths.MITOPORT_ALN, topology: Topology = paths.MITOPORT_TOPOLOGY,
    genes: GenesFile = paths.MITOPORT_GENES,
    species: Annotated[Path, typer.Option(help="the species of the matrix, one per line")] = paths.MITOPORT_SPECIES,
    threads: Annotated[int, typer.Option(help="MAFFT --thread; 4 as in the published run, 1 for a repeatable alignment")] = 4,
) -> None:
    """Align the nuclear proteins (MAFFT L-INS-i, from PATH or $MAFFT), filter, complete matrix, fixed topology."""
    from .mitoport import align
    align.main(data, sel, tree, outgroup, mafft_dir, aln_dir, topology, genes, species, threads)


@mitoport_app.command("branch-lengths")
def mitoport_branch_lengths(
    aln_dir: AlnDir = paths.MITOPORT_ALN, topology: Topology = paths.MITOPORT_TOPOLOGY,
    genes: GenesFile = paths.MITOPORT_GENES, bl_dir: BlDir = paths.MITOPORT_BL,
    threads: Annotated[int, typer.Option(help="IQ-TREE -T")] = 2,
    seed: Annotated[Optional[int], typer.Option(help="IQ-TREE --seed; a random one when omitted, as in the published run")] = None,
) -> None:
    """Branch lengths per protein on the fixed topology (IQ-TREE from PATH: -te, LG+G4 or mtMAM+G4, --keep-ident)."""
    from .mitoport import branch_lengths
    if branch_lengths.main(aln_dir, topology, genes, bl_dir, threads, seed):
        raise typer.Exit(1)


@mitoport_app.command("erc")
def mitoport_erc(
    genes: GenesFile = paths.MITOPORT_GENES, bl_dir: BlDir = paths.MITOPORT_BL,
    out: Annotated[Path, typer.Option("--out", help="the table to write")] = paths.MITOPORT_OUT / "erc.tsv",
) -> None:
    """Rate covariation between nuclear subunits and mitochondrial genes, against control proteins."""
    from .mitoport import erc
    erc.main(genes, bl_dir, out)


@mitoport_app.command("coembed")
def mitoport_coembed(
    embeddings: Annotated[Path, typer.Option(help="folder with <set>-esm2_650m.index.parquet and .mean.npy, from `download prepared --topic embeddings`")] = paths.EMBEDDINGS,
    tree: Annotated[Path, typer.Option(help="the dated species tree, from `longevity dated-tree`")] = paths.DATED_TREE,
    out: Out = paths.MITOPORT_OUT,
) -> None:
    """Do partners' ESM-2 embedding distances between species follow each other beyond the dated tree?"""
    from .mitoport import coembed
    coembed.main(embeddings, tree, out)


def main() -> None:
    """The console script: a missing input or program is one line and exit status 1, not a traceback."""
    from .external import MissingProgram
    try:
        app()
    except (FileNotFoundError, MissingProgram) as e:
        typer.echo(f"mitoage: {e}", err=True)
        raise SystemExit(1)
