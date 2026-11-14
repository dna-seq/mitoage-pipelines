# mitoage-pipelines

Analysis code behind the research arms of [MitoAge](https://just-dna.life/mitoage/): what the alignments of
the mitochondrial genes of about 2,300 animal species say about longevity, the peptides encoded inside
mtDNA, and the pairing of mitochondrial proteins with their nuclear partners (mitoport).

It is written so that a reader without access to our machines can rerun it: every input is downloaded by a
command of this repository, every step is one command, and every table a finding is read from is written
under `data/output`.

- **The pipeline in pictures**: [docs/PIPELINE.md](docs/PIPELINE.md).
- **What we asked, how we tested it and what came out**: [docs/HYPOTHESES.md](docs/HYPOTHESES.md). Results
  that did not hold, and one mistake of ours, are kept there next to those that did.
- The website and the pipeline that builds its data are in
  [dna-seq/mitoage](https://github.com/dna-seq/mitoage) (TypeScript). This repository is the Python side and
  reads the published build, [longevity-genie/mitoage](https://huggingface.co/datasets/longevity-genie/mitoage)
  on Hugging Face.

- **MitoAge build 2 from build 1**: [docs/BUILD2_REPRODUCTION.md](docs/BUILD2_REPRODUCTION.md). Every source
  of build 2 has a downloader here, every rule a step (`mitoage build2 ...`), and `mitoage build2 compare`
  checks the result against the published build 2: on the pinned inputs every value is identical (see
  "Build 2 reproduction" below).

Nothing here is peer reviewed.

## Setup

```sh
uv sync                                   # Python and its packages (uv.lock pins them)
conda env create -f environment.yml       # the external programs; or install them yourself
```

Python packages come from `uv`. The external programs do not: MAFFT 7.526, IQ-TREE 3.1.4, HyPhy 2.5.104
(and Coevol 1.6, built from its source). They are looked up on `PATH`; a step that needs one says which when
it is missing. `environment.yml` names the versions we ran; it has not been solved on a clean machine.

## Data layout

```
data/input     downloads, never edited: the MitoAge build, its dated species tree, NCBI ortholog sets,
               the sources of build 2 (data/input/build2);
               MANIFEST.json records what was fetched, from where, at which revision and when
data/interim   made on the way: alignments, trees, selections, inputs of other programs
data/output    the result tables
data/reference our own result tables, when fetched with `download prepared`, to compare a rerun with
```

All of it is ignored by git. NCBI's ortholog sets change with every annotation release, so the prepared
copy is the way to get the records we used. `MITOAGE_PIPELINES_DATA` moves the root elsewhere.

## Run

`uv run mitoage --help` lists everything; each command has `--help`. In order:

```sh
# inputs
uv run mitoage download build             # the MitoAge build from Hugging Face (about 46 MB)
uv run mitoage download tree              # the build's dated species tree, from the site's repository
uv run mitoage download orthologs         # nuclear subunits and control proteins from NCBI (62 genes)

# or: our prepared copies from Hugging Face, by topic (partners, embeddings, aligners, results, all)
uv run mitoage download prepared          # exactly what we fetched and computed; our tables land in data/reference

# H1: amino-acid content against lifespan
uv run mitoage longevity controls         # composition.tsv, and the earlier tree-free tests
uv run mitoage longevity fourfold
uv run mitoage longevity composition      # the earlier version (mixed model), kept as a lesson
uv run mitoage longevity dated-tree
uv run mitoage longevity aa-pgls --class Mammalia    # on the dated tree; also Aves, Actinopterygii, Reptilia

# H2: single residues against lifespan
uv run mitoage longevity sites            # the tree-free version, kept as a warning: it is inflated
uv run mitoage longevity pgls --class Mammalia
uv run mitoage longevity pgls --class Aves
uv run mitoage longevity pgls-power --class Mammalia
# on a second alignment (after `download prepared --topic aligners`); tools: famsa, muscle5, mafft_iter
uv run mitoage longevity pgls --class Mammalia --alignments-dir data/interim/aligners --tool famsa \
    --out data/output/longevity/pgls.aligner-famsa.Mammalia.tsv

# H3: selection in long-lived against short-lived lines (HyPhy)
uv run mitoage longevity rates
mitoage longevity motif-pgls --structure <mtdna_structure.parquet>
mitoage longevity relax-prep        # needs `longevity dated-tree`, like every step that takes a tree
uv run mitoage longevity relax

# H4: peptides encoded inside mtDNA
uv run mitoage peptides pilot

# H5: mitochondrial and nuclear partners (MAFFT, IQ-TREE)
uv run mitoage mitoport select
uv run mitoage mitoport align
uv run mitoage mitoport branch-lengths
uv run mitoage mitoport erc

# H6: embeddings of the partners (needs `download prepared --topic embeddings`)
uv run mitoage mitoport coembed
```

The longest steps are RELAX (5 to 40 minutes per gene on 4 cores) and the tree of a class (minutes); the
rest take seconds to a few minutes.

## Build 2 reproduction

MitoAge build 2 (2026-11-13, 2,335 species) was made by the site's TypeScript pipeline from build 1 (the 2015
data, 922 species) and current upstream sources. `mitoage build2` makes it again from the same sources, with
the rules written as Python, and compares the result with the published tables. What build 1 and build 2
are, every source with its pin and licence, every rule in plain words with the site file that defines it,
the run times and the comparison are in [docs/BUILD2_REPRODUCTION.md](docs/BUILD2_REPRODUCTION.md).

```sh
uv run mitoage download hf                          # build 1 (build1/) and the published build 2, revision 653afe2a3b
uv run mitoage download refseq --from-snapshot      # RefSeq release 237: NCBI keeps only the current release
uv run mitoage download taxonomy --from-snapshot    # the NCBI taxonomy dump build 2 read (2026-11-13)
uv run mitoage download lifespan                    # AnAge, FishBase/SeaLifeBase v26.06, USGS BBL; tries EURING
uv run mitoage download cited && uv run mitoage download contig && uv run mitoage download legacy
uv run mitoage build2 refseq && uv run mitoage build2 match && uv run mitoage build2 status
uv run mitoage build2 primary-list && uv run mitoage download primary
uv run mitoage build2 assemble                      # species, region_stats, codon_usage, genome_features,
                                                    # lifespan_observations, changelog, review_list
uv run mitoage download tree-sources --from-snapshot && uv run mitoage build2 tree
uv run mitoage build2 align                         # MAFFT v7.526 from PATH
uv run mitoage download tree --ref 8ff350c504bb59c9171eaed729eb6cddaada0f0b
uv run mitoage build2 compare                       # -> data/output/build2/compare
```

About seven minutes in all. On 2026-11-14, with these inputs, all 5,302,404 values of the seven tables are
identical to the published build 2, the 28 alignment files have the same SHA-256 and the dated tree is
byte-identical. With today's sources instead (RefSeq release 238, today's taxonomy, live Open Tree of Life
and TimeTree) the species set is the same and the only change of values is one record NCBI re-annotated in
place (*Pongo pygmaeus* NC_001646.1); the document lists every difference with its reason.

Where it comes from:

- The pipeline that made build 2: [dna-seq/mitoage](https://github.com/dna-seq/mitoage) — its
  [README](https://github.com/dna-seq/mitoage/blob/main/README.md),
  [docs/BUILD2_REPORT.md](https://github.com/dna-seq/mitoage/blob/main/docs/BUILD2_REPORT.md),
  [docs/DECISIONS.md](https://github.com/dna-seq/mitoage/blob/main/docs/DECISIONS.md) and
  [plans/DATA_BUILD2_DECISIONS.md](https://github.com/dna-seq/mitoage/blob/main/plans/DATA_BUILD2_DECISIONS.md).
- The build: [longevity-genie/mitoage](https://huggingface.co/datasets/longevity-genie/mitoage) (build 2 at
  revision `653afe2a3b79e18ae892e579ded9acdb17e0c334`, build 1 in `build1/`).
- The inputs no upstream archive keeps (RefSeq release 237, the taxonomy dump of 2026-11-13, the GenBank
  records read by efetch, AnAge build 15, the USGS table, the Open Tree and TimeTree answers, the old site's
  pages): [longevity-genie/mitoage-build2-sources](https://huggingface.co/datasets/longevity-genie/mitoage-build2-sources),
  revision `9d55d72648d0838c3b5f37fca6185f8b938b3272`. Every downloader reads the original source by
  default and this snapshot with `--from-snapshot`; no token is needed. `mitoage snapshot build|upload` is
  for maintainers only (it needs a write token in `HF_TOKEN`).

## How far it reproduces

Checked on 2026-11-14 by running the commands above against the outputs of the original runs:

- **Identical bytes**: `controls`, `fourfold`, `composition`, `sites`, `pgls` (mammals and birds, on the
  dated tree and on a gene-built tree), `peptides pilot`, `relax-prep`, `coevol-prep`, `mitoport select`.
- **Equal to the second decimal**: `mitoport erc` (compared on the gene-built tree's topology, as first run). MAFFT with several threads and IQ-TREE without a seed are
  not repeatable run to run; the group medians move by at most 0.002. `--threads 1` and `--seed` make them so.
- **Does not reproduce**: HyPhy RELAX. On identical inputs a second run of ATP8 gave K = 1.90, p = 0.20
  where the first gave K = 1.49, p = 0.55: the optimiser lands on different optima. Read every RELAX number
  in the hypotheses document with that in mind; it is one more reason we claim nothing from that test.
- **Not rerun**: `download orthologs` for all 62 genes (one gene was downloaded and matched the earlier
  download record for record); `relax` for 12 of the 13 genes; the aligner comparison for 11 of 13 genes.

## Known gaps

- A tree built from the 13 proteins themselves (`longevity concat`, `longevity tree`, then `--gene-tree`,
  or `--tree` with `--outgroup`) is kept only to show how much a result depends on the tree. Every step
  takes the build's dated species tree by default.
- The embeddings (ESM-2 650M) are downloaded, not computed here: the script that made them is in another
  project. Only the mean over residues is published; the per-residue embeddings (about 8 GB) are not.
- The aligner comparison (`longevity compare-aligners`) compares alignment files you give it; the commands
  that produced them with each aligner are not wrapped.
- Lifespan as a continuous trait (Coevol): the inputs are written, the run and its reading are not wrapped.

## Tests

```sh
uv run pytest
```

Small real cases with expected values computed by hand in the test: the FASTA reader, fourfold-site
counts, the tree covariance of the regression, H1's regression against Felsenstein's contrasts on five
mammals and its species filter, H2's power formula; and the rules of build 2 (`tests/test_build2.py`: the base-label
fix of build 1, join() segments, the lifespan precedence and its tie-break, name matching, id minting, the
change classes, record_status against the changelog), each checked to fail when the code it guards is broken.

## Licence

Code: MIT. The data it downloads carries the licences named on the dataset's card and by NCBI.
