"""The `mitoage` command: every step, with the folders of `paths.py` as defaults.

Order of a full run:
  mitoage download build && mitoage download orthologs
  mitoage longevity controls && mitoage longevity fourfold && mitoage longevity composition
  mitoage longevity sites
  mitoage download tree && mitoage longevity dated-tree && mitoage longevity aa-pgls && mitoage longevity pgls && mitoage longevity pgls-power
  mitoage longevity rates && mitoage longevity motif-pgls
  mitoage longevity relax-prep && mitoage longevity relax
  mitoage longevity coevol-prep
  mitoage peptides pilot
  mitoage mitoport select && mitoage mitoport align && mitoage mitoport branch-lengths && mitoage mitoport erc
  mitoage download prepared --topic embeddings && mitoage mitoport coembed
  mitoage download pdb 5Z62 && mitoage mitoport contacts && mitoage mitoport burial
  mitoage download prepared --topic residues && mitoage mitoport iface-align && mitoage mitoport iface
  mitoage download pdb 5XTD 5Z62 8H9S 5XTE 5XTH 5XTI 6ZQN 7AJD && mitoage structure within
  mitoage structure fold (GPU), or mitoage download prepared --topic structures; then mitoage structure msa
  mitoage longevity coevol --chain 1 && mitoage longevity coevol --chain 2 && mitoage longevity coevol-read
Build 2 from build 1 (docs/BUILD2_REPRODUCTION.md):
  mitoage download hf|refseq|taxonomy|lifespan|cited|contig|legacy (--from-snapshot where no archive keeps it)
  mitoage build2 refseq && mitoage build2 match && mitoage build2 status
  mitoage build2 primary-list && mitoage download primary && mitoage build2 assemble
  mitoage download tree-sources && mitoage build2 tree && mitoage build2 align && mitoage build2 compare
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
structure_app = typer.Typer(no_args_is_help=True, help="Structural alignments of the 13 mitochondrial proteins: a check of the sequence alignments.")
snapshot_app = typer.Typer(no_args_is_help=True, help="MAINTAINERS ONLY: build and upload the snapshot dataset of build 2's sources (needs a write token).")
build2_app = typer.Typer(no_args_is_help=True, help="MitoAge build 2 from build 1 and the upstream sources (docs/BUILD2_REPRODUCTION.md).")
app.add_typer(download_app, name="download")
app.add_typer(longevity_app, name="longevity")
app.add_typer(peptides_app, name="peptides")
app.add_typer(mitoport_app, name="mitoport")
app.add_typer(structure_app, name="structure")
app.add_typer(build2_app, name="build2")
app.add_typer(snapshot_app, name="snapshot")

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
    topic: Annotated[str, typer.Option(help="partners, embeddings, residues (7.7 GB, not in all), aligners, structures, results, or all")] = "all",
    revision: Annotated[Optional[str], typer.Option(help="a commit of the dataset, to pin what is read (one topic only)")] = None,
) -> None:
    """What we fetched and computed, from Hugging Face, by topic: start any step without redoing those before it.
    Our result tables go to data/reference, for comparison with a rerun's data/output."""
    from . import download
    topics = [t for t in download.PREPARED if t != "residues"] if topic == "all" else [topic]
    if topic != "all" and topic not in download.PREPARED:
        raise typer.BadParameter(f"unknown topic {topic!r}; one of {', '.join(download.PREPARED)} or all")
    for t in topics:
        download.prepared(t, paths.ROOT, paths.MANIFEST, revision if len(topics) == 1 else None)


@download_app.command("tree")
def download_tree(ref: Annotated[str, typer.Option(help="branch, tag or commit of the site's repository")] = "main") -> None:
    """The build's dated species tree (phylogeny.json) from the site's repository into data/input/mitoage-site."""
    from . import download
    download.site_tree(paths.SITE_TREE, paths.MANIFEST, ref)


@download_app.command("pdb")
def download_pdb(ids: Annotated[list[str], typer.Argument(help="PDB entry ids, e.g. 5Z62")]) -> None:
    """Solved structures as mmCIF from RCSB into data/input/pdb, with each file's SHA-256 in the manifest."""
    from . import download
    download.pdb([i.upper() for i in ids], paths.PDB, paths.MANIFEST)


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


# --- download: the sources of build 2 -----------------------------------------------------------
Manifest = Annotated[Path, typer.Option(help="the download record")]


Snapshot = Annotated[bool, typer.Option("--from-snapshot", help="read the copy in the snapshot dataset longevity-genie/mitoage-build2-sources (no token) instead of the original source")]
SnapRev = Annotated[Optional[str], typer.Option("--snapshot-revision", help="revision of the snapshot dataset; the pinned one when omitted")]


def _snapshot(kind: str, revision: Optional[str], manifest: Path) -> None:
    from .build2 import pins, snapshot
    snapshot.download(kind, revision or pins.SNAPSHOT_REVISION, manifest)


@download_app.command("hf")
def download_hf(revision: Annotated[str, typer.Option(help="commit of the dataset; build 2's when omitted")] = "653afe2a3b79e18ae892e579ded9acdb17e0c334",
                dest: Annotated[Path, typer.Option(help="where it goes")] = paths.B2_HF, manifest: Manifest = paths.MANIFEST) -> None:
    """MitoAge build 1 (build1/) and the published build 2 tables and alignments, at a pinned revision."""
    from .build2 import download
    download.hf(dest, manifest, revision)


@download_app.command("refseq")
def download_refseq(
    expect: Annotated[Optional[str], typer.Option(help="refuse unless NCBI's current release is this one")] = None,
    from_file: Annotated[Optional[Path], typer.Option("--from-file", help="a kept copy of release 237, accepted only with build 2's md5")] = None,
    efetch: Annotated[bool, typer.Option("--efetch", help="pinned mode without the release file: exactly build 2's accession.versions, by E-utilities")] = False,
    from_snapshot: Snapshot = False, snapshot_revision: SnapRev = None,
    manifest: Manifest = paths.MANIFEST,
) -> None:
    """The RefSeq mitochondrion release (NCBI keeps only the current one), release 237 from the snapshot, a kept copy, or --efetch."""
    from .build2 import download
    if from_snapshot:
        _snapshot("refseq", snapshot_revision, manifest)
    elif from_file:
        download.refseq_from_file(from_file, paths.B2_REFSEQ, manifest)
    elif efetch:
        import polars as pl
        sp = paths.B2_HF / "species.parquet"
        if not sp.is_file():
            raise FileNotFoundError(f"no {sp}: run `mitoage download hf` first (the version list is build 2's refseq column)")
        versions = pl.read_parquet(sp).filter(pl.col("record_status") != "retired")["refseq"].to_list()
        download.refseq_efetch(versions, paths.B2_REFSEQ, manifest, label="build2-653afe2a3b")
    else:
        download.refseq_current(paths.B2_REFSEQ, manifest, expect)


@download_app.command("taxonomy")
def download_taxonomy(
    archive: Annotated[Optional[str], typer.Option(help="date of a monthly archive (YYYY-MM-01); 2026-11-01 is the closest to build 2's dump of 2026-11-13")] = None,
    current: Annotated[bool, typer.Option("--current", help="the current new_taxdump instead of an archive")] = False,
    from_file: Annotated[Optional[Path], typer.Option("--from-file", help="a kept new_taxdump.tar.gz of 2026-11-13, accepted only with build 2's md5")] = None,
    from_snapshot: Snapshot = False, snapshot_revision: SnapRev = None,
    manifest: Manifest = paths.MANIFEST,
) -> None:
    """NCBI Taxonomy (new_taxdump): a monthly archive, the current dump, build 2's dump of 2026-11-13 from the snapshot, or a kept copy."""
    from .build2 import download
    if from_snapshot:
        _snapshot("taxonomy", snapshot_revision, manifest)
    elif from_file:
        download.taxdump_from_file(from_file, paths.B2_TAXDUMP, manifest)
    elif current:
        download.taxdump_current(paths.B2_TAXDUMP, manifest)
    else:
        download.taxdump_archive(archive or "2026-11-01", paths.B2_TAXDUMP, manifest)


@download_app.command("lifespan")
def download_lifespan(manifest: Manifest = paths.MANIFEST,
                      fishbase_version: Annotated[str, typer.Option(help="FishBase/SeaLifeBase snapshot on source.coop")] = "v26.06",
                      from_snapshot: Snapshot = False, snapshot_revision: SnapRev = None) -> None:
    """AnAge, FishBase and SeaLifeBase (source.coop), USGS Bird Banding Lab longevity records; tries EURING."""
    from .build2 import download
    if from_snapshot:
        _snapshot("lifespan", snapshot_revision, manifest)
        download.euring(paths.B2_EURING, manifest)
        return
    download.anage(paths.B2_ANAGE, manifest)
    download.fishbase("fb", paths.B2_FISHBASE, manifest, fishbase_version)
    download.fishbase("slb", paths.B2_SEALIFEBASE, manifest, fishbase_version)
    download.usgs(paths.B2_USGS, manifest)
    download.euring(paths.B2_EURING, manifest)


@download_app.command("cited")
def download_cited(manifest: Manifest = paths.MANIFEST, from_snapshot: Snapshot = False, snapshot_revision: SnapRev = None) -> None:
    """The RefSeq records build 1 cites (one accession.version per species), by E-utilities efetch."""
    if from_snapshot:
        return _snapshot("cited", snapshot_revision, manifest)
    from .build2 import build1, download
    raw = build1.load(paths.B2_BUILD1)
    download.genbank([s.refseq for s in raw.kept], paths.B2_GENBANK / "cited", manifest, "the records cited by build 1")


@download_app.command("legacy")
def download_legacy(manifest: Manifest = paths.MANIFEST, from_snapshot: Snapshot = False, snapshot_revision: SnapRev = None) -> None:
    """The old site's species ids, from its total_mtDNA export pages (www.mitoage.info)."""
    if from_snapshot:
        return _snapshot("legacy", snapshot_revision, manifest)
    from .build2 import download
    download.legacy(paths.B2_LEGACY, manifest)


@download_app.command("contig")
def download_contig(release: Annotated[str, typer.Option(help="the release whose CONTIG-only records are fetched: a number, or efetch")] = "237",
                    manifest: Manifest = paths.MANIFEST, from_snapshot: Snapshot = False, snapshot_revision: SnapRev = None) -> None:
    """Sequence of the complete metazoan records that carry only a CONTIG line in the release, by efetch."""
    if from_snapshot:
        return _snapshot("contig", snapshot_revision, manifest)
    from .build2 import download, run
    download.contig(run.release_file(release)[0], paths.B2_GENBANK / "contig", manifest)


@download_app.command("primary")
def download_primary(manifest: Manifest = paths.MANIFEST,
                     species: Annotated[Path, typer.Option(help="hagrids whose cited record's INSDC source is fetched, one per line; written by `build2 primary-list`")] = paths.B2_INTERIM / "primary_hagrids.txt",
                     from_snapshot: Snapshot = False, snapshot_revision: SnapRev = None) -> None:
    """The INSDC records named in the COMMENT of build 1's cited RefSeq records (for the re-annotation test)."""
    if from_snapshot:
        return _snapshot("primary", snapshot_revision, manifest)
    from .build2 import download
    download.primary(species, paths.B2_GENBANK / "cited", paths.B2_GENBANK / "primary", manifest)


@download_app.command("tree-sources")
def download_tree_sources(species: Annotated[Path, typer.Option(help="species.parquet of the build the tree is for")] = paths.B2_OUT / "species.parquet",
                          manifest: Manifest = paths.MANIFEST, from_snapshot: Snapshot = False, snapshot_revision: SnapRev = None) -> None:
    """The Open Tree of Life and TimeTree answers the dated tree needs (asks only for those not yet in data/input/build2/tree).
    With --from-snapshot: the answers build 2's tree was built from."""
    if from_snapshot:
        return _snapshot("tree", snapshot_revision, manifest)
    from .build2 import download, tree
    t = tree.main(species, paths.B2_TREE_IN, paths.B2_INTERIM / "phylogeny.download-check.json", online=True)
    download._manifest(manifest, "opentreeoflife+timetree", f"{tree.pins.OTL_API} induced_subtree, tnrs; {tree.pins.TIMETREE_API} mrca",
                       paths.B2_TREE_IN, revision=f"OTL synthetic tree {t['version']}; TimeTree public API (unversioned)")


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


@longevity_app.command("rates")
def longevity_rates(
    data: Build = paths.BUILD, composition: Composition = COMPOSITION, fourfold: Fourfold = FOURFOLD,
    tree: Annotated[Path, typer.Option("--tree", help="the dated species tree, from `longevity dated-tree`")] = paths.DATED_TREE,
    out: Out = paths.LONGEVITY_OUT,
) -> None:
    """H7: substitution rates against lifespan, by the relative-rate test on sister pairs of the dated tree."""
    from .longevity import rates
    rates.main(data, composition, fourfold, tree, out)


@longevity_app.command("motif-pgls")
def longevity_motif_pgls(
    data: Build = paths.BUILD, composition: Composition = COMPOSITION, fourfold: Fourfold = FOURFOLD,
    tree: Annotated[Path, typer.Option("--tree", help="the dated species tree, from `longevity dated-tree`")] = paths.DATED_TREE,
    out: Out = paths.LONGEVITY_OUT,
    structure: Annotated[Optional[Path], typer.Option("--structure", help="mtdna_structure.parquet from another build, while the motif features are not yet in a released one")] = None,
) -> None:
    """H8: the build's mtDNA motif, repeat and contact-zone features against lifespan (H1's model)."""
    from .longevity import motif_pgls
    motif_pgls.main(data, composition, fourfold, tree, out, structure)


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


@longevity_app.command("aa-pgls")
def longevity_aa_pgls(
    data: Build = paths.BUILD, composition: Composition = COMPOSITION, fourfold: Fourfold = FOURFOLD,
    cls: Class = "Mammalia",
    tree: Annotated[Optional[Path], typer.Option(help="the tree; the dated species tree when omitted")] = None,
    out: Annotated[Optional[Path], typer.Option("--out", help="the table to write; <output>/longevity/aa_pgls.<Class>.tsv when omitted")] = None,
) -> None:
    """H1 on the dated tree: content of each amino acid against lifespan, mass and fourfold composition (Pagel's lambda)."""
    from .longevity import aa_pgls
    aa_pgls.main(data, composition, fourfold, tree or paths.DATED_TREE, out or paths.aa_pgls_table(cls), cls)


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


@longevity_app.command("pgls-power")
def longevity_pgls_power(cls: Class = "Mammalia",
                         table: Annotated[Optional[Path], typer.Option(help="the column test's table; pgls.<Class>.tsv when omitted")] = None,
                         out: Annotated[Optional[Path], typer.Option("--out", help="the table to write; pgls_power.<Class>.tsv when omitted")] = None) -> None:
    """H2's power: the smallest residue effect on lifespan detectable at 80 % power and the alpha 5 % FDR needs."""
    from .longevity import power
    power.main(table or paths.pgls_table(cls), out or paths.pgls_power_table(cls))


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


@longevity_app.command("coevol")
def longevity_coevol(
    chain: Annotated[int, typer.Option(help="chain number; the published run had chains 1 and 2")] = 1,
    cycles: Annotated[int, typer.Option(help="cycles to save; 1,100 in the published run")] = 1100,
    prep: Annotated[Path, typer.Option(help="folder of `coevol-prep`'s inputs")] = paths.COEVOL,
    out: Out = paths.COEVOL_OUT,
) -> None:
    """One Coevol chain (dS and dN/dS jointly with lifespan and mass; coevol from PATH or $COEVOL). Days of CPU."""
    from .longevity import coevol
    coevol.run(prep, out, chain, cycles)


@longevity_app.command("coevol-read")
def longevity_coevol_read(
    chains: Annotated[list[int], typer.Option("--chain", help="chains to read (repeat the option)")] = [1, 2],
    burnin: Annotated[int, typer.Option(help="cycles discarded; 200 in the published reading")] = 200,
    out: Out = paths.COEVOL_OUT,
) -> None:
    """readcoevol per chain and tracecomp across chains: correlations, and whether the chains have converged."""
    from .longevity import coevol
    coevol.read(out, chains, burnin)


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


CifOpt = Annotated[Path, typer.Option("--cif", help="the solved structure of complex IV, mmCIF")]
IfaceDir = Annotated[Path, typer.Option("--work", help="folder of the contact map: contacts.tsv, sites.tsv, chains.tsv, species.txt")]


@mitoport_app.command("contacts")
def mitoport_contacts(cif: CifOpt = paths.PDB / "5Z62.cif", data: Build = paths.BUILD, out: IfaceDir = paths.MITOPORT_IFACE) -> None:
    """Residues of COX1-3 within 4.5 A of another subunit in the solved complex IV (after `download pdb 5Z62`)."""
    from .mitoport import contacts
    contacts.main(cif, data, out)


@mitoport_app.command("burial")
def mitoport_burial(
    cif: CifOpt = paths.PDB / "5Z62.cif", data: Build = paths.BUILD, work: IfaceDir = paths.MITOPORT_IFACE,
    tree: Annotated[Path, typer.Option(help="the dated species tree, from `longevity dated-tree`")] = paths.DATED_TREE,
    quality: Annotated[Path, typer.Option(help="species with a usable lifespan record, from `longevity dated-tree`")] = paths.QUALITY_OK,
) -> None:
    """Accessibility of each COX1-3 residue, the burial-matched non-contact control, columns and occupancy."""
    from .mitoport import burial
    burial.main(cif, work / "contacts.tsv", data, tree, quality, work)


@mitoport_app.command("iface-align")
def mitoport_iface_align(
    sel: Sel = paths.MITOPORT_SEL, out: Annotated[Path, typer.Option("--out", help="folder for <GENE>.afa")] = paths.MITOPORT_IFACE_ALN,
    threads: Annotated[int, typer.Option(help="MAFFT --thread per gene; 2 as published, 1 for a repeatable alignment")] = 2,
    jobs: Annotated[int, typer.Option(help="genes aligned at once")] = 8,
) -> None:
    """Every nuclear selection with MAFFT L-INS-i (--anysymbol), unfiltered: the alignments the interface test reads."""
    from .mitoport import iface
    iface.align(sel, out, threads, jobs)


@mitoport_app.command("iface")
def mitoport_iface(
    cif: CifOpt = paths.PDB / "5Z62.cif", data: Build = paths.BUILD, work: IfaceDir = paths.MITOPORT_IFACE,
    nucaln: Annotated[Path, typer.Option("--nucaln", help="the alignments of `iface-align`")] = paths.MITOPORT_IFACE_ALN,
    esm: Annotated[Path, typer.Option(help="folder with {mt13,nuclear}-esm2_650m.index.parquet and .residues.npy, from `download prepared --topic residues`")] = paths.EMBEDDINGS,
    skip_esm: Annotated[bool, typer.Option("--skip-esm", help="run without the embeddings: every esm_ column is left empty (not assessed)")] = False,
    tree: Annotated[Path, typer.Option(help="the dated species tree, from `longevity dated-tree`")] = paths.DATED_TREE,
    out: Out = paths.IFACE_OUT,
) -> None:
    """Coupling of contact residues with their nuclear partners over species pairs, and lifespan on the contact mismatch."""
    from .mitoport import iface
    iface.main(work, nucaln, work / "contacts.tsv", cif, data, esm, tree, out, skip_esm)


# --- structure --------------------------------------------------------------------------------
@structure_app.command("within")
def structure_within(
    candidates: Annotated[Optional[list[str]], typer.Argument(help="PDB entries tried in order; the published list when omitted")] = None,
    data: Build = paths.BUILD, out: Out = paths.STRUCTURE_OUT,
) -> None:
    """The solved human proteins of the complexes superimposed all against all (US-align from PATH)."""
    from .structure import within
    within.main(data, paths.PDB, candidates or within.CANDIDATES, paths.STRUCTURE_CHAINS, out)


@structure_app.command("fold")
def structure_fold(
    data: Build = paths.BUILD,
    species: Annotated[Path, typer.Option(help="species to fold, one per line; the 62 mammals of `mitoport align`")] = paths.MITOPORT_SPECIES,
    out: Annotated[Path, typer.Option("--out", help="folder for <GENE>/<species>.pdb and plddt.tsv")] = paths.STRUCTURE_PDBS,
    checkpoint: Annotated[str, typer.Option(help="ESMFold2 checkpoint on Hugging Face")] = "biohub/ESMFold2-Fast",
    samples: Annotated[int, typer.Option(help="diffusion samples per protein; 1 as published")] = 1,
) -> None:
    """Predicted structures of the 13 proteins with ESMFold2 (needs a CUDA GPU and `uv sync --group fold`)."""
    from .structure import fold
    fold.main(data, species, out, checkpoint, None, samples)


@structure_app.command("msa")
def structure_msa(
    pdbs: Annotated[Path, typer.Option(help="the predicted structures of `structure fold`")] = paths.STRUCTURE_PDBS,
    data: Build = paths.BUILD, out: Out = paths.STRUCTURE_OUT,
    threads: Annotated[int, typer.Option(help="MAFFT --thread")] = 8,
) -> None:
    """FoldMason's structural alignment against MAFFT L-INS-i on the same sequences: agreement and cost."""
    from .structure import msa
    msa.main(pdbs, data, paths.STRUCTURE_TMP, out, None, threads)


# --- build 2 ----------------------------------------------------------------------------------
Release = Annotated[str, typer.Option(help="RefSeq input: 237 (build 2's release), another release number, or efetch (build 2's versions by E-utilities)")]
Taxonomy = Annotated[str, typer.Option(help="taxonomy dump folder under data/input/build2/taxdump: current-2026-11-13 (build 2's), 2026-11-01 (archive), current-<date>")]
B2Out = Annotated[Path, typer.Option("--out", help="folder the tables are written to")]


@build2_app.command("refseq")
def build2_refseq(release: Release = "237", taxonomy: Taxonomy = "current-2026-11-13") -> None:
    """Complete metazoan RefSeq mitogenomes, one per NCBI species (refseq_records.parquet)."""
    from .build2 import run
    run.stage_refseq(run.inputs(release, taxonomy), paths.B2_INTERIM)


@build2_app.command("match")
def build2_match(release: Release = "237", taxonomy: Taxonomy = "current-2026-11-13") -> None:
    """Lifespan-source names matched to NCBI species (observations_matched.parquet)."""
    from .build2 import run
    run.stage_match(run.inputs(release, taxonomy), paths.B2_INTERIM)


@build2_app.command("status")
def build2_status(release: Release = "237", taxonomy: Taxonomy = "current-2026-11-13") -> None:
    """Status of the 922 build 1 species and the candidate new species (build1_status.parquet, candidates.parquet)."""
    from .build2 import run
    run.stage_status(run.inputs(release, taxonomy), paths.B2_INTERIM)


@build2_app.command("primary-list")
def build2_primary_list(release: Release = "237", taxonomy: Taxonomy = "current-2026-11-13") -> None:
    """The build 1 species whose INSDC source record `download primary` fetches (primary_hagrids.txt)."""
    from .build2 import run
    run.stage_primary_list(run.inputs(release, taxonomy), paths.B2_INTERIM)


@build2_app.command("assemble")
def build2_assemble(release: Release = "237", taxonomy: Taxonomy = "current-2026-11-13", out: B2Out = paths.B2_OUT) -> None:
    """The build 2 tables: species, region_stats, codon_usage, genome_features, lifespan_observations, changelog, review_list."""
    from .build2 import run
    run.stage_assemble(run.inputs(release, taxonomy), out)


@build2_app.command("tree")
def build2_tree(species: Annotated[Path, typer.Option(help="species.parquet of the build")] = paths.B2_OUT / "species.parquet",
                out: Annotated[Path, typer.Option("--out", help="the tree to write")] = paths.B2_OUT / "phylogeny.json") -> None:
    """The dated species tree (Open Tree of Life topology, TimeTree ages) from the answers in data/input/build2/tree; offline."""
    from .build2 import tree
    tree.main(species, paths.B2_TREE_IN, out, online=False)


@build2_app.command("align")
def build2_align(release: Release = "237", out: B2Out = paths.B2_OUT / "alignments",
                 species: Annotated[Path, typer.Option(help="species.parquet of the build whose live species are aligned")] = paths.B2_OUT / "species.parquet",
                 threads: Annotated[int, typer.Option(help="MAFFT --thread; 4 as in the published run")] = 4) -> None:
    """Alignments of the 13 proteins (and their codons) and the two rRNAs with MAFFT v7.526 FFT-NS-2 (from PATH or $MAFFT)."""
    from .build2 import align, run
    rf, rl = run.release_file(release)
    m = align.main(species, rf, rl, paths.B2_GENBANK / "contig", out, threads)
    print(f"{len(m['files'])} alignment files for {m['species']} species -> {out}")


@build2_app.command("compare")
def build2_compare(reproduced: Annotated[Path, typer.Option(help="the reproduction's tables")] = paths.B2_OUT,
                   published: Annotated[Path, typer.Option(help="the published build 2 (download hf)")] = paths.B2_HF,
                   site_tree: Annotated[Path, typer.Option(help="build 2's dated tree, from the site repository (download tree)")] = paths.SITE_TREE,
                   out: Annotated[Optional[Path], typer.Option("--out", help="where the report goes; <reproduced>/compare when omitted")] = None) -> None:
    """Agreement of the reproduction with the published build 2, per table and field; every disagreement with a reason."""
    from .build2 import compare
    compare.main(published, reproduced, out or reproduced / "compare", site_tree=site_tree)


# --- snapshot dataset (maintainers) ---------------------------------------------------------------
Staging = Annotated[Path, typer.Option(help="the folder the snapshot is assembled in")]


@snapshot_app.command("build")
def snapshot_build(
    release237: Annotated[Path, typer.Option(help="RefSeq release 237 as build 2 read it (md5-checked)")] = paths.refseq_release("237"),
    taxdump_tar: Annotated[Path, typer.Option(help="new_taxdump.tar.gz of 2026-11-13 as build 2 read it (md5-checked)")] = Path("new_taxdump.tar.gz"),
    tree_answers: Annotated[Path, typer.Option(help="folder of the OTL/TimeTree answers build 2's tree read")] = paths.B2_TREE_IN,
    tree_used: Annotated[Path, typer.Option(help="phylogeny_answers.txt written by `build2 tree`: the answers it read")] = paths.B2_OUT / "phylogeny_answers.txt",
    dates: Annotated[str, typer.Option(help="JSON: retrieval dates (UTC) of release237, efetch, anage, fishbase, usgs, tree, legacy")] = "{}",
    staging: Staging = paths.B2_INTERIM / "snapshot",
) -> None:
    """Assemble the snapshot dataset (files, MANIFEST.json, card) from data/input."""
    import json
    from .build2 import snapshot
    d = json.loads(dates)
    used = [l for l in tree_used.read_text().split() if l]
    e = snapshot.build(staging, paths.B2_IN, {"release237": release237, "taxdump_tar": taxdump_tar, "tree_answers": tree_answers, "tree_used": used,
                                              **{f"{k}_retrieved": v for k, v in d.items()}})
    print(f"{len(e)} files, {sum(x['bytes'] for x in e):,} bytes in {staging}")


@snapshot_app.command("upload")
def snapshot_upload(staging: Staging = paths.B2_INTERIM / "snapshot") -> None:
    """Upload the snapshot as one commit (HF_TOKEN must be set in this command's environment; it is never printed)."""
    from .build2 import snapshot
    print(f"{snapshot.REPO}: commit {snapshot.upload(staging)}")


def main() -> None:
    """The console script: a missing input or program is one line and exit status 1, not a traceback."""
    from .external import MissingProgram
    from .build2.download import Refused as RefusedDownload
    try:
        app()
    except (FileNotFoundError, MissingProgram, RefusedDownload) as e:
        typer.echo(f"mitoage: {e}", err=True)
        raise SystemExit(1)
