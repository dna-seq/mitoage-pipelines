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
  on Hugging Face. Every dataset this repository downloads is listed in the
  [MitoAge collection](https://huggingface.co/collections/longevity-genie/mitoage-6af952ac7256661d479577a5) there.

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

Python packages come from `uv`. The external programs do not: MAFFT 7.526, IQ-TREE 3.1.4, HyPhy 2.5.104,
US-align 20241108 and FoldMason 4.dd3c235 (and Coevol 1.6, built from its source at commit 2661e95 of
github.com/bayesiancook/coevol). Folding (`structure fold`) needs a CUDA GPU and the `fold` dependency group,
which brings esm and torch (the published run: esm 3.4.3, torch 2.15.0, Python 3.14). They are looked up on `PATH`; a step that needs one says which when
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

# or: our prepared copies from Hugging Face, by topic (partners, embeddings, residues, aligners, structures,
# results, all; `all` leaves out the 7.7 GB of per-residue embeddings, `--topic residues` fetches them)
# The coevolution data is there for reuse beyond our tests: `partners` holds the NCBI ortholog sets and alignments
# of the nuclear subunits (mammals; complex IV across vertebrates); `results` holds rate covariation (erc.tsv),
# embedding covariation (coembed.tsv), complex IV's contact map and residue coupling (output/mitoport/iface/), and
# every species' mito-nuclear pair discordance (partners_species.tsv).
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
uv run mitoage download technology && uv run mitoage longevity technology   # H8's control: each record's sequencing technology, from NCBI
uv run mitoage longevity motif-pgls      # H8: the build's motif, repeat and contact-zone features, four classes
uv run mitoage longevity structure-pgls  # H12: complex IV contact mismatch (after `mitoport burial`) and GRAVY, four classes
# H12, second round: the partners of complex IV in four classes, and ESM-2 at its contacts (GPU, `uv sync --group fold`)
uv run mitoage download orthologs --ortholog-taxon vertebrata --dest data/input/ncbi-vertebrata --gene COX4I1 ...  # the 11 subunits
uv run mitoage mitoport select --ncbi data/input/ncbi-vertebrata --out data/interim/mitoport/partners/sel \
    --class Mammalia --class Aves --class Actinopterygii --class Reptilia
uv run mitoage mitoport iface-align --sel data/interim/mitoport/partners/sel --out data/interim/mitoport/partners/aln --threads 1
# or: `download prepared --topic partners` brings those ortholog sets, selections and alignments as published
uv run mitoage mitoport partners          # discordance of the mito-nuclear contact pairs against lifespan
uv run mitoage structure embed            # ESM-2 650M means at the contact and control residues (a CUDA GPU)
uv run mitoage longevity structure-pgls --family embedding --out data/output/longevity/embedding
uv run mitoage longevity generation-time   # H11: lifespan, maturity, gestation, birth weight and the slow-fast axis; mammals and birds (on haumea: slurm/generation-time.sbatch)
mitoage longevity relax-prep        # needs `longevity dated-tree`, like every step that takes a tree
uv run mitoage longevity relax

# H4: peptides encoded inside mtDNA
uv run mitoage peptides pilot
# H15 (RM16): the whole catalogue across the build, with the window, synonymous and dinucleotide nulls, the tree and
# lifespan (after `longevity dated-tree`, `controls`, `fourfold`; on haumea: slurm/job4-peptides.sbatch)
uv run mitoage peptides catalogue && uv run mitoage peptides frames
uv run mitoage peptides nulls --jobs 8 && uv run mitoage peptides tree && uv run mitoage peptides pgls

# H5: mitochondrial and nuclear partners (MAFFT, IQ-TREE)
uv run mitoage mitoport select
uv run mitoage mitoport align
uv run mitoage mitoport branch-lengths
uv run mitoage mitoport erc

# H6: embeddings of the partners (needs `download prepared --topic embeddings`)
uv run mitoage mitoport coembed

# H5 at the residue level and H9: the contacts of complex IV (after `mitoport select` and `longevity dated-tree`)
uv run mitoage download pdb 5Z62                  # the human enzyme from RCSB; its SHA-256 goes to the manifest
uv run mitoage mitoport contacts                  # residue pairs within 4.5 A
uv run mitoage mitoport burial                    # accessibility, the burial-matched control, columns
uv run mitoage mitoport iface-align               # every nuclear selection, MAFFT L-INS-i
uv run mitoage download prepared --topic residues   # the per-residue embeddings, 7.7 GB
uv run mitoage mitoport iface                     # or --skip-esm: the esm_ columns left empty, not assessed

# H14: the contact maps of complexes I, III and V (H5's design and H12's pairs; on haumea: slurm/job1-contacts.sbatch)
uv run mitoage download pdb 5XTD 5XTE 8H9S
uv run mitoage mitoport contacts --complex I && uv run mitoage mitoport burial --complex I   # also III and V
uv run mitoage download orthologs --gene NDUFV1 ...          # complex I and V subunits (genes.COMPLEX_I_NUCLEAR, COMPLEX_V_NUCLEAR)
uv run mitoage mitoport iface --complex I --skip-esm          # mammals; no per-residue embeddings for the new genes
uv run mitoage mitoport partners --complex I                  # vertebrate partners of the subunits that touch
uv run mitoage mitoport readouts                              # the readouts fixed in advance, every complex

# H19: the coupling class by class, against distance from the mtDNA-encoded core, and its mitochondrial side
# (after the contact maps of the four complexes; on haumea: slurm/h19-classes.sbatch)
uv run mitoage download orthologs --ortholog-taxon vertebrata --dest data/input/ncbi-vertebrata --gene NDUFV1 ...   # genes.H19_NUCLEAR
uv run mitoage mitoport select --ncbi data/input/ncbi-vertebrata --out data/interim/mitoport/partners/sel \
  --class Mammalia --class Aves --class Actinopterygii --class Reptilia
uv run mitoage mitoport iface-align --sel data/interim/mitoport/partners/sel --out data/interim/mitoport/partners/aln --threads 1
uv run mitoage mitoport classes                               # core_distance, class_coupling, class_matrix, class_readouts

# H16: the published cybrids against a mismatch count at the mito-nuclear contacts (seconds on a laptop)
uv run mitoage download pdb 5XTD 5XTE 5Z62 8H9S
uv run mitoage download cybrid                    # 15 RefSeq mitogenomes at fixed versions and 3 cell-line cytb records
uv run mitoage mitoport cybrid                    # positions, scores, calibration and summary tables

# H17: the cybrids with scores that weigh substitutions (on haumea: slurm/h17-cybrid.sbatch, then h17-cybrid-esm.sbatch)
uv run mitoage mitoport contacts --complex I && uv run mitoage mitoport burial --complex I   # and III, IV, V (as H14)
uv run mitoage mitoport cybrid-coev               # co-evolved residue pairs the chimera loses, contact against control
uv sync --group fold && uv run mitoage mitoport cybrid-esm   # ESM-2 masked marginals, alone and with the partner (GPU)
uv run mitoage mitoport cybrid-readout            # H16's readouts for these scores, beside H16's baselines

# H18: the cybrids by an interface energy on the chimeric complexes (on haumea: slurm/h18-energy.sbatch, CPU nodes)
uv sync --group energy                            # OpenMM and PDBFixer
uv run mitoage download pdb 5XTD 5XTE 5Z62 8H9S 8OM1 7O3H 7O37
uv run mitoage mitoport cybrid-energy --list      # the 44 tasks; one with --task N, or --host, --donor and --complex
uv run mitoage mitoport cybrid-energy-collect && uv run mitoage mitoport cybrid-energy-readout

# H10: structural alignments
uv run mitoage download pdb 5XTD 5Z62 8H9S 5XTE 5XTH 5XTI 6ZQN 7AJD
uv run mitoage structure within                   # US-align, the solved human proteins all against all
uv sync --group fold && uv run mitoage structure fold   # ESMFold2 on a CUDA GPU: 62 mammals x 13 proteins
# or: uv run mitoage download prepared --topic structures   # the 806 structures of the published run
uv run mitoage structure msa                      # FoldMason against MAFFT on the same sequences

# H13: an audit of the build's records (strand, gene order, genetic code, duplicates and mislabels; not a longevity test)
uv run mitoage download refseq --from-snapshot && uv run mitoage download contig --from-snapshot
uv run mitoage download taxonomy --from-snapshot  # RefSeq release 237 and the taxonomy build 2 read
uv run mitoage audit records                      # every species' record from the pinned release; inputs.tsv: every input's SHA-256
uv run mitoage audit strand                       # conventional, reversed, mixed or unknown, beside the GC/AT skews
uv run mitoage audit order                        # gene order canonicalised by rotation; arrangements per class, order, family; flags
uv run mitoage audit code                         # declared against NCBI Taxonomy's genetic code; internal stops under each
uv run mitoage audit duplicates                   # pairs at >= 99 % identity; nearest neighbours outside their order or class
# (on haumea: slurm/job3-audit.sbatch, which also runs the tests, a mutation check and a rerun)

# H3's continuous-trait model: Coevol 1.6, built from source (days of CPU per chain)
uv run mitoage longevity coevol --chain 1 && uv run mitoage longevity coevol --chain 2
uv run mitoage longevity coevol-read              # readcoevol and tracecomp; read only once the chains converge
uv run mitoage longevity coevol-resume --chain 1 --stop-after-hours 70   # continue a chain, checkpoint every 6 h
uv run mitoage longevity coevol-prep --fixed-ages && uv run mitoage longevity coevol --fixed-ages --chain 1
                                                  # the same with the dated tree's ages held (-fixtimes)
uv run mitoage longevity branch-rates             # the two-step check: dS and dN/dS per branch (HyPhy, hours)
uv run mitoage longevity branch-pgls              # then terminal branches against lifespan and mass
```

On haumea, every analysis the MitoAge manuscript reports runs on one build as one Slurm job,
`slurm/manuscript-rerun.sbatch` (its header gives the commands); `slurm/compare-runs.py` compares two such runs
table by table and class by class.

The longest steps are RELAX (5 to 40 minutes per gene on 4 cores) and the tree of a class (minutes); the
rest take seconds to a few minutes.

## Build 2 reproduction

MitoAge build 2 (2026-11-13, with the EURING bird records added on 2026-11-15: 2,344 species) was made by the site's TypeScript pipeline from build 1 (the 2015
data, 922 species) and current upstream sources. `mitoage build2` makes it again from the same sources, with
the rules written as Python, and compares the result with the published tables. What build 1 and build 2
are, every source with its pin and licence, every rule in plain words with the site file that defines it,
the run times and the comparison are in [docs/BUILD2_REPRODUCTION.md](docs/BUILD2_REPRODUCTION.md).

```sh
uv run mitoage download hf                          # build 1 (build1/) and the published build 2, revision 0f4e970554
uv run mitoage download refseq --from-snapshot      # RefSeq release 237: NCBI keeps only the current release
uv run mitoage download taxonomy --from-snapshot    # the NCBI taxonomy dump build 2 read (2026-11-13)
uv run mitoage download lifespan                    # AnAge, FishBase/SeaLifeBase v26.06, USGS BBL, EURING (Internet Archive)
uv run mitoage download cited && uv run mitoage download contig && uv run mitoage download legacy
uv run mitoage build2 refseq && uv run mitoage build2 match && uv run mitoage build2 status
uv run mitoage build2 primary-list && uv run mitoage download primary
uv run mitoage build2 assemble                      # species, region_stats, codon_usage, genome_features,
                                                    # lifespan_observations, changelog, review_list
uv run mitoage download tree-sources --from-snapshot && uv run mitoage build2 tree
uv run mitoage build2 align                         # MAFFT v7.526 from PATH
uv run mitoage download tree --ref 03c09498564c17d1445bc4c18c70d0e710d471db
uv run mitoage build2 compare                       # -> data/output/build2/compare
```

About six minutes in all. On 2026-11-15, with these inputs, all 5,328,887 values of the seven tables are
identical to the published build 2, the 28 alignment files have the same SHA-256 and the dated tree is
byte-identical. On the build before EURING (2,335 species, 2026-11-14), with today's sources instead (RefSeq
release 238, today's taxonomy, live Open Tree of Life and TimeTree) the species set was the same and the only change of values is one record NCBI re-annotated in
place (*Pongo pygmaeus* NC_001646.1); the document lists every difference with its reason.

Where it comes from:

- The pipeline that made build 2: [dna-seq/mitoage](https://github.com/dna-seq/mitoage) — its
  [README](https://github.com/dna-seq/mitoage/blob/main/README.md),
  [docs/BUILD2_REPORT.md](https://github.com/dna-seq/mitoage/blob/main/docs/BUILD2_REPORT.md),
  [docs/DECISIONS.md](https://github.com/dna-seq/mitoage/blob/main/docs/DECISIONS.md) and
  [plans/DATA_BUILD2_DECISIONS.md](https://github.com/dna-seq/mitoage/blob/main/plans/DATA_BUILD2_DECISIONS.md).
- The build: [longevity-genie/mitoage](https://huggingface.co/datasets/longevity-genie/mitoage) (build 2 at
  revision `0f4e970554b2fe0265e718321f422347293dde17`, the build with the EURING birds and the hand-checked gene
  annotations, build 1 in `build1/`).
- The inputs no upstream archive keeps (RefSeq release 237, the taxonomy dump of 2026-11-13, the GenBank
  records read by efetch, AnAge build 15, the USGS table, the Open Tree and TimeTree answers, the old site's
  pages): [longevity-genie/mitoage-build2-sources](https://huggingface.co/datasets/longevity-genie/mitoage-build2-sources),
  revision `a76215144ec33d6ed1e5ea599cb797f2651dab1a`. The EURING PDF is not in it: the Internet Archive keeps it.
  Every downloader reads the original source by default and this snapshot with `--from-snapshot`; no token is
  needed. `mitoage snapshot build|upload` is for maintainers only (a write token in `HF_TOKEN`, or the machine's
  stored Hugging Face login).

## How far it reproduces

Checked on 2026-11-14 by running the commands above against the outputs of the original runs:

- **Identical bytes**: `controls`, `fourfold`, `composition`, `sites`, `pgls` (mammals and birds, on the
  dated tree and on a gene-built tree), `peptides pilot`, `relax-prep`, `coevol-prep`, `mitoport select`.
- **Identical bytes on a rerun of 2026-11-15**: `longevity generation-time` (both maturity modes; Slurm job 6487335 on haumea, from 858974b; tables in `docs/analyses/h11/`).
- **Identical bytes on a rerun, on haumea, 2026-11-15**: `longevity structure-pgls` (H12); `longevity motif-pgls`
  gives the table H8 reports byte for byte.
- **Identical bytes on a rerun, on haumea, 2026-11-15 (H12's second round)**: `longevity structure-pgls` (with the
  own-composition model, and `--family embedding` on one embedding table) and `mitoport partners`. `structure embed`
  itself (ESM-2 on a GPU) was run once.
- **Identical bytes on a rerun, on haumea, 2026-11-15, of everything the manuscript reports**
  (`slurm/manuscript-rerun.sbatch`, job 6489876 from 2ee27a1, on the build before the EURING birds): `longevity
  rates` (`rates.tsv`, `rate_contrasts.tsv`), `motif-pgls`, `generation-time` in both modes, and the rows of
  `structure-pgls` that `results/h12/` holds. The same job on the build with the EURING birds (job 6489877) gave the
  tables in `results/euring/`, whose README lists what the nine birds and the realignment changed.
- **Identical bytes on a rerun, on haumea, 2026-11-15**: every `audit` step (H13; job 6490112 from 49d52d5, build
  c525adb; tables and log in `results/job3-audit/`; and job 6490297 from fb1cf00 on build 0f4e970, the manuscript's,
  in `results/job3-audit/0f4e970/`).
- **Equal to the second decimal**: `mitoport erc` (compared on the gene-built tree's topology, as first run). MAFFT with several threads and IQ-TREE without a seed are
  not repeatable run to run; the group medians move by at most 0.002. `--threads 1` and `--seed` make them so.
- **Does not reproduce**: HyPhy RELAX. On identical inputs a second run of ATP8 gave K = 1.90, p = 0.20
  where the first gave K = 1.49, p = 0.55: the optimiser lands on different optima. Read every RELAX number
  in the hypotheses document with that in mind; it is one more reason we claim nothing from that test.
- **The steps first run on the cluster on 2026-11-14 and 15** (by another agent, as one-off scripts, then
  ported here) were rerun from this code on 2026-11-15 on the same cluster and compared with that run's files:
  byte-identical for `longevity dated-tree`, `download pdb 5Z62`, `mitoport contacts`, `mitoport burial`
  (sites, chains, species), `mitoport select` (62 of 62 selections) and `longevity coevol-prep`; the same
  tables for `mitoport iface` (coupling.tsv identical on that run's nuclear alignments, lifespan.tsv identical
  with the interval columns added), `structure within` and `structure msa` (on that run's structures, timings
  aside). `mitoport iface-align` gave 60 of 62 alignments identical: MAFFT with two threads is not repeatable,
  and on the fresh alignments one protein's sequence correlations move (UQCRQ, +0.565 to +0.572), which
  changes no group median; `--threads 1` makes the alignments repeatable.
- **Not rerun**: `structure fold` (3.3 GPU-hours) and the Coevol chains (days of CPU each); the `esm_` columns
  of `mitoport iface` were rerun on the cluster with the per-residue embeddings, which are now on Hugging Face
  (`--topic residues`; their SHA-256 is on the dataset's card), and `structure msa` on the structures now
  published as `--topic structures`.
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
- Lifespan as a continuous trait (Coevol): `longevity coevol` and `coevol-read` wrap the run and its reading;
  the published chains did not converge (effective sample size 7 for dN/dS), so its numbers are provisional.
- The per-residue ESM-2 embeddings that `mitoport iface` reads for its `esm_` columns (7.7 GB) are downloaded,
  not computed here: like the mean embeddings, they were made in another project. Without them `--skip-esm`
  leaves those columns empty, marked not assessed; the `seq_` columns and the lifespan test are unaffected.
- The 806 predicted structures of H10 (150 MB) are downloaded with `--topic structures` or made again by
  `structure fold` on a GPU (about 3.3 GPU-hours on an A100); a second folding is not expected to give the
  same coordinates to the last digit, since the model samples.

## Tests

```sh
uv run pytest
```

Small real cases with expected values computed by hand in the test: the FASTA reader, the complex IV
contact map's pieces (subunit names, column mapping, the burial-matched control, the distances and the
modal mismatch, the regression), the structural check's pair agreement and core, the US-align and Coevol
readers on real output kept in `tests/data`, fourfold-site
counts, the tree covariance of the regression, H1's regression against Felsenstein's contrasts on five
mammals and its species filter, H2's power formula; and the rules of build 2 (`tests/test_build2.py`: the base-label
fix of build 1, join() segments, the lifespan precedence and its tie-break, name matching, id minting, the
change classes, record_status against the changelog), each checked to fail when the code it guards is broken; and the
audit of H13 (`tests/test_audit.py`, on real records of RefSeq release 237: a reverse-complemented record swaps
conventional and reversed, a rotated one keeps its arrangement, edited bovids are flagged with the kind of edit, the
translation under the declared code is RefSeq's own, identical records pair; `slurm/job3-audit-mutations.py` breaks
25 guards one at a time).
`tests/test_complexes.py` (H14) builds the contact maps of complexes III, IV and V from the structures, checks that
complex IV's is byte for byte the published one, and tests the readouts; `slurm/job1-contacts-mutations.py` breaks
each of its guards in turn.
`tests/test_cybrid.py` (H16) checks the gene names, the placing of a protein on human positions (a partial
sequence's ends are unknown, not deleted), the counts and the ROC area by hand, the agreement of sources, the shipped
measurement table, and that the contact residues of all four complexes are the set `mitoport contacts` finds; each
test was checked to fail when the code it guards was broken (the cutoff, the second copy of cytochrome b, the
partial ends, the disagreement of sources, the radical substitutions).

`tests/test_cybrid_h17.py` (H17) works the sequence weights, the joint distribution, the mutual information and the loss
of a pair out by hand, checks that host and donor are left out of their own counts and that unknown or thin pairs are
unscored and never zero, and tests the ESM-2 bookkeeping on a table of log-probabilities in place of the model (which
positions are masked, pairs over the model's length, the partner's effect); `slurm/h17-cybrid-mutations.py` breaks 11
guards one at a time.

`tests/test_cybrid_energy.py` (H18) checks by hand which chain is which protein, what is substituted and what is
counted as not modelled, the tasks and the collection; with OpenMM installed, the cutting of a chain at a break, the
free and restrained atoms, and on COX2 of 5Z62 the system, the substitutions made and that two groups moved apart do
not interact; `slurm/h18-energy-mutations.py` breaks 12 guards one at a time.

## Licence

Code: MIT. The data it downloads carries the licences named on the dataset's card and by NCBI.
