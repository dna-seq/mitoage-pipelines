# H37 part A: the 22 tRNA and 2 rRNA genes against maximum lifespan

Run of 2026-11-26. The bars and nulls are those of entry H37 in `docs/HYPOTHESES.md` as committed in 0d65412, unchanged. Part B (where in the
tRNA) and part C (the tRNA against its own synthetase) are not run. The code is commit f9f8fdd; the tables below are those of two runs from
fresh clones of it that are identical byte for byte.

## In one paragraph

**Part A does not count: under a correction over all 768 cells every cell is null, and no tRNA or rRNA gene is named.** The
paragraph below is the narrower reading first run, under the entry's own correction within each test; the wider one, which
supersedes it, is under "The correction over all 768 cells".

768 cells were judged (24 genes, each in 7 composition tests and 1 rate test in each of 4 taxon groups). Under the entry's
own correction (Holm over the 24 genes within one test) **2 count**, both
of one gene in one group: *tRNA-Thr* of ray-finned fishes, its GC skew and its AT skew, which are two readings of one difference in how much T and C the
gene carries. 4 composition cells and 3 rate cells pass the first condition alone and are **not counting**, in the entry's words:
the GC content, the C share and the T share of the 16S rRNA in mammals and the A share of tRNA-Gln in mammals, and the rate of tRNA-Cys, tRNA-Ser(AGY)
and tRNA-Tyr in mammals; each fails the third condition, the high-confidence subset. The other 666 composition cells and 93 rate
cells are null, none of them uninformative (0), and none is not assessed (0). The positive control reproduced H7's table exactly, so the
rates are used. What the two counting cells do and do not show is under "Reading".

## Inputs

- The build: Hugging Face dataset `longevity-genie/mitoage` at 446dad61883937075ec0c7e1f3f8135f2962a1cf, the revision H27 ran on (it carries each
  species' confidence level, which H27's `high` set needs). Species table SHA-256 fa967e33724af62dab5859eab48393404627e54aee951bcd7fe69895b86a6db1.
  Dated species tree: `phylogeny.json` of `dna-seq/mitoage` at e214d8e1dbbaae836f5708f2d5efbac307e8c9d4, SHA-256
  bfcd086948b3e6d44b485201e1cd1af4c536360eebc935fd851d9d3c2ae19e76 (the same bytes as at 03c0949, the tree of H7's published numbers).
- The records: NCBI RefSeq mitochondrion release 237, md5 4491799651740f49c8f6a011cf9a96e5 (the pin of build 2, checked by the code), and the sequence of
  the nine CONTIG-only records from `download contig --from-snapshot`: 2325 live species, every one with its record.
- The aligner: MAFFT v7.526, `--retree 2 --maxiterate 0 --quiet`, the rule of the build's own nucleotide alignments (`build2/align.py`), one
  process and one thread per gene, species sorted by slug. The two rRNA alignments are the build's, not realigned: the sequences extracted from
  the records equal the build's rows for all 2323 (12S) and 2322 (16S) species (`rrna_check.tsv`).
- For the control on the build of H7's published numbers: the species table of build 0f4e970, SHA-256
  267a16311cca579af2ca1cd004e99befe66e53e1ae7696d71ffb9f6b7be5dbed, with the same tree. `composition.tsv` (SHA-256
  d8786a480349ffe6502192d8e09f4bd6fa49dd30d1f0d3ab437af1b8fd025f67) is the same bytes on both builds.

## The positive control, before any tRNA number was read

The rate driver of this part, in codon mode, on the 13 protein genes and their concatenation, with H7's species (the default dataset of H1, body
mass, a tip of the tree, all 13 proteins), cherries, outgroups, Nei-Gojobori counting and tests, was run on the build of H7's published numbers and on the
build of this run, and compared row by row, at the five significant digits the table is written at, with H7's `results/review-20261118/longevity/rates.tsv`
(the corrected run of 2026-11-18). **655 rows each, 0 rows and 0 cells differ, none only on one side.** Both written tables are byte-identical to
the published one (SHA-256 b2b656ad2ea5cfa17ec22579e364636e39433cee78a8535c6b360f28733e6791), and so is H27's table of its unrestricted run on this
build. Among the numbers the entry quotes: mammals, with the body-mass contrast in the model, concatenated genes: slope +0.112, p = 0.0222,
q = 0.155, 115 pairs; 340 cherries in all. `h7_control.0f4e970.json` and `h7_control.446dad6.json` hold the two verdicts. The code refuses
`rna-rates` and `rna-tests` without a reproduced control of the same build.

What this control covers: the cherries, the outgroups, the contrasts, the regression, the saturation flag and the correction, which the RNA rates share
with H7. It does not cover the distance itself: an RNA has no dN or dS, and the Jukes-Cantor distance on shared columns is new code, tested on real
tRNA sequences and by deliberate breaks, but not reproduced against any earlier number.

## What each species carries

Per gene, before any test (`gene_census.tsv`, `alignments.tsv`, `unresolved.tsv`). 2325 live species. A gene is the first copy in genome order, in its own
sense.

| gene | species with it | not annotated | unresolved isoacceptor | annotated twice | informative columns, alignment | informative columns, smallest group |
|---|---|---|---|---|---|---|
| tRNA-Ala | 2322 | 3 | 0 | 3 | 68 | 49 |
| tRNA-Arg | 2322 | 3 | 0 | 0 | 67 | 50 |
| tRNA-Asn | 2323 | 2 | 0 | 3 | 73 | 50 |
| tRNA-Asp | 2322 | 3 | 0 | 7 | 69 | 58 |
| tRNA-Cys | 2322 | 3 | 0 | 1 | 65 | 42 |
| tRNA-Gln | 2321 | 4 | 0 | 8 | 71 | 45 |
| tRNA-Glu | 2319 | 6 | 0 | 23 | 69 | 55 |
| tRNA-Gly | 2322 | 3 | 0 | 1 | 67 | 32 |
| tRNA-His | 2321 | 4 | 0 | 2 | 68 | 42 |
| tRNA-Ile | 2320 | 5 | 0 | 5 | 69 | 34 |
| tRNA-Leu(CUN) | 2279 | 4 | 42 | 10 | 72 | 29 |
| tRNA-Leu(UUR) | 2270 | 11 | 44 | 3 | 73 | 22 |
| tRNA-Lys | 2320 | 5 | 0 | 6 | 70 | 39 |
| tRNA-Met | 2323 | 2 | 0 | 22 | 66 | 24 |
| tRNA-Phe | 2319 | 6 | 0 | 8 | 69 | 49 |
| tRNA-Pro | 2314 | 11 | 0 | 26 | 69 | 44 |
| tRNA-Ser(AGY) | 2272 | 13 | 40 | 3 | 65 | 47 |
| tRNA-Ser(UCN) | 2277 | 8 | 40 | 8 | 71 | 33 |
| tRNA-Thr | 2319 | 6 | 0 | 26 | 69 | 28 |
| tRNA-Trp | 2323 | 2 | 0 | 1 | 69 | 25 |
| tRNA-Tyr | 2322 | 3 | 0 | 3 | 71 | 23 |
| tRNA-Val | 2319 | 6 | 0 | 3 | 66 | 50 |
| rRNA-12S | 2323 | 2 | 0 | 7 | 905 | 544 |
| rRNA-16S | 2322 | 3 | 0 | 7 | 1483 | 963 |

- *Not annotated*: the record has no feature of that name. *Unresolved*: the record has a Leu or Ser tRNA whose isoacceptor neither the annotation nor the
  vertebrate gene order (chordates only) places, and the gene is missing from that record: 94 such Leu copies and 78 Ser copies in all, left out of
  both alignments of their amino acid and counted here. *Annotated twice*: the first copy is used.
- *Informative columns*: columns with a base in at least half the rows and two different bases each in at least two rows, over the whole alignment and, in
  the last column, the smallest over the four taxon groups' species of the tests. The rRNAs have 905 and 1483 over the alignment.
- Twenty tRNA features are not 50 to 100 nucleotides long (nine shorter, eleven longer; the longest are a tRNA-Ile of 3539 nucleotides in
  *Chanodichthys mongolicus*, one of 1026 in *Otolemur crassicaudatus* and a tRNA-Gln of 611 in *Synthliboramphus antiquus*). They are annotations of the
  records, and are kept; they make the alignments of tRNA-Ile (3601 columns), tRNA-Gln (667) and tRNA-Asp (436) wider than a tRNA and do not change their
  informative columns.

## What was fixed before the run, where the entry is silent

These are in the code from commit 3028d9d, before any number existed.

- *The features*: GC content, GC skew and AT skew as the build's region statistics compute them (`genome.make_stats`), and the shares of A, C, G and T among
  the resolved bases (the analogue of H1's fourfold shares), each on the gene in its own sense. H1's amino-acid and fourfold features need codons, which an RNA
  has not; the rate is the other measure.
- *The rate*: H7's cherry contrast with the Jukes-Cantor distance on the columns A, C, G or T in all three species (H7's shared-columns rule), at least 30 such
  columns (H7's floor is 30 codons); saturated when a distance exceeds 1 (H7's flag); the three normalisations that apply (the raw difference, per million years
  of the cherry's depth, with the body-mass contrast). The bar reads the one with the body mass in the model, saturated pairs left out, H7's headline. The
  protein genes' counterpart is H7's dS with the body-mass contrast.
- *The composition test*: per taxon group (Mammalia, Aves, Actinopterygii, Reptilia), a phylogenetic regression on the dated tree of the standardised
  feature on log10 lifespan and log10 adult body mass (`aa_pgls.fit`, H1's and H11's), with Pagel's lambda estimated (the primary model) and with lambda fixed at
  1, the model under which the regression equals the one on independent contrasts; at least 25 species. H1's species. The covariance is the group's, a gene using the
  rows of its own species.
- *The high-confidence subset*: H27's `high` (confidence level high), applied to the same species before the regression or the cherries.
- *The effect*: the absolute partial correlation of lifespan, computed from the t statistic and its degrees of freedom, the same for a regression and for a
  contrast; a measure that does not depend on the units of the gene.
- *The bar, as coded*: (1) the Holm-adjusted p over the 24 genes, within each set, group and test, is below 0.05 with lambda estimated and with lambda 1, with the
  same sign (a rate is a contrast on cherries, the contrasts control by construction); (2) the effect is at least the median of the 13 protein genes' in the same
  group, test and set; (3) on the high set, Holm-adjusted p below 0.05 with the same sign. Holm, because the entry says family-wise; Benjamini-Hochberg's q is
  in the tables beside it. The median needs at least 7 of the 13 protein genes.
- *Few informative columns*: fewer than 20 among the species of the test; the entry gives no number.

## Results

768 cells: **2 count, 4 composition and 3 rate cells are not counting, 666 composition and 93 rate cells
are null; 0 are uninformative nulls and 0 are not assessed.** The smallest number of informative columns in any cell is 22. Every gene,
group and test has its row in `rna_verdict.tsv`; `rna_gene_summary.tsv` has one row per gene.

The nine cells that pass the first condition. Slopes are in standard deviations of the feature per tenfold lifespan for a composition test, and a raw contrast
slope for a rate; the effect is the absolute partial correlation; "protein median" is the median of the 13 protein genes' in the same group and test.

| test | group | gene | n | informative columns | slope | Holm p, lambda estimated | Holm p, lambda 1 | effect (absolute r) | protein median | n, high | slope, high | Holm p, high | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| composition: A_share | Mammalia | tRNA-Gln | 603 | 60 | -0.747 | 0.0367 | 0.0364 | 0.129 | 0.052 | 527 | -0.712 | 0.115 | not counting |
| composition: C_share | Mammalia | rRNA-16S | 603 | 1098 | +0.474 | 0.00381 | 0.00379 | 0.153 | 0.081 | 527 | +0.392 | 0.121 | not counting |
| composition: T_share | Mammalia | rRNA-16S | 603 | 1098 | -0.464 | 0.0264 | 0.0264 | 0.133 | 0.092 | 527 | -0.392 | 0.321 | not counting |
| composition: at_skew | Actinopterygii | tRNA-Thr | 235 | 58 | -0.372 | 0.0148 | 0.0034 | 0.222 | 0.063 | 184 | -0.387 | 0.0235 | counts |
| composition: gc_percent | Mammalia | rRNA-16S | 603 | 1098 | +0.485 | 0.00395 | 0.00393 | 0.153 | 0.098 | 527 | +0.409 | 0.0956 | not counting |
| composition: gc_skew | Actinopterygii | tRNA-Thr | 235 | 58 | +0.469 | 0.0122 | 6.04e-05 | 0.226 | 0.102 | 184 | +0.541 | 0.00647 | counts |
| rate: rate | Mammalia | tRNA-Cys | 165 | 59 | +0.178 | 0.0142 | 0.0142 | 0.264 | 0.083 | 123 | +0.128 | 0.511 | not counting |
| rate: rate | Mammalia | tRNA-Ser(AGY) | 163 | 56 | -0.071 | 0.031 | 0.031 | 0.248 | 0.083 | 123 | -0.071 | 0.0722 | not counting |
| rate: rate | Mammalia | tRNA-Tyr | 167 | 57 | +0.154 | 0.00965 | 0.00965 | 0.271 | 0.083 | 125 | +0.105 | 0.424 | not counting |

Reading the table: in ray-finned fishes the GC skew of tRNA-Thr rises and its AT skew falls with lifespan, at an effect 2.2 times (GC skew) and 3.5 times (AT skew) the protein-gene
median, with both lambda models and on the high set (Holm p 0.0065 and 0.023), and **counts**. In mammals the GC content, the C share and the T share of the 16S rRNA, the
A share of tRNA-Gln and the rate of tRNA-Cys, tRNA-Ser(AGY) and tRNA-Tyr pass the controls and are larger than the protein median, and each is **not counting**: on the
high-confidence subset the same sign remains but the Holm-adjusted p is 0.072 to 0.51. The mammalian 16S rRNA GC content agrees in direction with H11's GC content of the whole record (rising with the slow
life history); that is a resemblance, not tested here.

The medians of the 13 protein genes' effects, which the bar compares with (default species set, lambda estimated; for the rate, dS with the body-mass contrast):

| test | Mammalia | Aves | Actinopterygii | Reptilia |
|---|---|---|---|---|
| A_share | 0.052 | 0.038 | 0.080 | 0.055 |
| C_share | 0.081 | 0.025 | 0.036 | 0.060 |
| G_share | 0.082 | 0.037 | 0.081 | 0.045 |
| T_share | 0.092 | 0.028 | 0.036 | 0.061 |
| at_skew | 0.061 | 0.042 | 0.063 | 0.047 |
| gc_percent | 0.098 | 0.024 | 0.064 | 0.072 |
| gc_skew | 0.025 | 0.030 | 0.102 | 0.065 |
| rate | 0.083 | 0.057 | 0.085 | not computed |

The reptiles' rate has no protein-gene median (no single protein gene has 12 sister pairs with dS at or below 1 in the reptiles; the unsaturated pairs are the bar's), so a reptile tRNA that had passed the first condition would
have been "not assessed (no protein-gene median)"; none did.

Per gene, over its 32 cells: cells passing (1), cells that count, its smallest Holm p, and its largest effect over the protein median (1 is the median) with the cell it
is in. The largest ratio is a maximum over 32 cells and is large for most genes by chance alone; it is a place against the median, not a result.

| gene | cells passing (1) | cells that count | smallest Holm p | largest effect over the protein median | in the cell |
|---|---|---|---|---|---|
| tRNA-Ala | 0 | 0 | 0.189 | 5.7 | composition at_skew Reptilia |
| tRNA-Arg | 0 | 0 | 1 | 3.8 | composition at_skew Reptilia |
| tRNA-Asn | 0 | 0 | 1 | 3.9 | composition gc_percent Aves |
| tRNA-Asp | 0 | 0 | 0.277 | 3.4 | composition C_share Actinopterygii |
| tRNA-Cys | 1 | 0 | 0.0142 | 5.0 | composition T_share Actinopterygii |
| tRNA-Gln | 1 | 0 | 0.0367 | 2.6 | composition C_share Actinopterygii |
| tRNA-Glu | 0 | 0 | 0.156 | 4.3 | composition gc_skew Reptilia |
| tRNA-Gly | 0 | 0 | 1 | 3.0 | composition C_share Reptilia |
| tRNA-His | 0 | 0 | 0.548 | 4.5 | composition G_share Reptilia |
| tRNA-Ile | 0 | 0 | 0.557 | 4.8 | composition G_share Reptilia |
| tRNA-Leu(CUN) | 0 | 0 | 1 | 2.3 | composition G_share Reptilia |
| tRNA-Leu(UUR) | 0 | 0 | 0.958 | 3.9 | rate rate Actinopterygii |
| tRNA-Lys | 0 | 0 | 0.146 | 5.9 | composition at_skew Reptilia |
| tRNA-Met | 0 | 0 | 0.271 | 3.3 | composition G_share Reptilia |
| tRNA-Phe | 0 | 0 | 0.594 | 2.9 | composition gc_percent Aves |
| tRNA-Pro | 0 | 0 | 0.88 | 3.9 | composition gc_percent Aves |
| tRNA-Ser(AGY) | 1 | 0 | 0.031 | 3.0 | rate rate Mammalia |
| tRNA-Ser(UCN) | 0 | 0 | 0.509 | 3.2 | composition G_share Aves |
| tRNA-Thr | 2 | 2 | 0.0122 | 4.2 | composition gc_percent Aves |
| tRNA-Trp | 0 | 0 | 0.191 | 4.1 | composition gc_percent Aves |
| tRNA-Tyr | 1 | 0 | 0.00965 | 3.3 | rate rate Mammalia |
| tRNA-Val | 0 | 0 | 0.124 | 4.4 | rate rate Actinopterygii |
| rRNA-12S | 0 | 0 | 0.376 | 3.3 | composition G_share Reptilia |
| rRNA-16S | 3 | 0 | 0.00381 | 4.2 | composition gc_skew Mammalia |

## The correction over all 768 cells

Added 2026-11-26, code 8ddbb66, on the same tables (no job). The entry's correction is within each test, over the 24 genes:
it holds the chance of a false pass in *one* test's family at 5%. There are 32 such families (7 composition features and the
rate in 4 groups), so a chance pass in at least one of them under no effect at all has probability 1 - 0.95^32, about 0.81,
were they independent. The same bar was therefore read again against Holm over the whole family of 768 cells of a species
set at once, each lambda model on its own (`correct_all`, `judge_wide`; the `_all` columns and `verdict_all` of
`rna_verdict.tsv`).

**Every one of the 768 cells is null under it: 0 count, 0 are not counting.** What the two cells that counted become:

| cell (tRNA-Thr, ray-finned fishes) | Holm within its test | Holm over all 768 | BH q over all 768 |
| --- | --- | --- | --- |
| GC skew, lambda estimated | 0.012 | 0.387 | 0.079 |
| GC skew, lambda 1 | 0.00006 | 0.0019 | 0.0019 |
| AT skew, lambda estimated | 0.015 | 0.471 | 0.079 |
| AT skew, lambda 1 | 0.0034 | 0.108 | 0.017 |

The bar's first condition needs both lambda models and the primary one, lambda estimated, fails in both cells. The second
condition (the effect against the protein-gene median) is the same number as before, corrected for nothing, and the third
is read against the wider correction too.

Three things say this is not an artefact of the family's exact size:

- With lambda estimated, the smallest p of all 768 cells is 1.6e-4 (the mammalian 16S rRNA's C share), Holm 0.122; **no cell
  anywhere reaches q below 0.05 under Benjamini-Hochberg either**, which is the lenient correction. Under lambda 1 alone 14
  cells do, 13 of them tRNA-Thr's own and its neighbours' in the same group.
- A defensible count of independent cells cannot be given exactly: the 7 composition features of a gene are functions of its
  four base shares, which sum to 1, so a gene and group hold at most 3 free composition dimensions, and the genes of a group
  share one tree and one set of species. At the most generous reduction that bound allows, 24 genes x 3 x 4 groups plus the
  96 rates, 384 cells, tRNA-Thr's GC skew is 0.19 with lambda estimated and the best cell of the whole family is 0.061:
  still nothing passes.
- The two cells are one signal (one difference in T and C content), so counting them as two never helped the result.

## Reading

- **Part A names no gene.** Under the entry's own correction one gene in one group passes, tRNA-Thr in ray-finned fishes, in
  two cells that are one signal; under the correction over all 768 cells it does not, and the wider reading is the one that
  holds. The bar was not changed to fit either result: both are reported, and the entry says the result does not count.
- Not examined: whether tRNA-Thr's annotation in fishes is uniform (26 records of the build annotate it twice); whether a few clades of fishes carry the
  signal; any mechanism. These are unknown.
- **No tRNA or rRNA gene shows a rate that follows lifespan** under the bar: the three mammalian tRNA rates that pass the first condition fail the third.
- **How small an effect these tests could have seen is unknown**: no power calculation was made. A tRNA gene has 22 to 68 informative columns (median 53) among the
  species of a group, so the nulls are weak in the entry's sense even where none is labelled uninformative.

## Checks

- Tests: 64 of 64 of this part pass with none skipped; breaks: 87 of 87 caught (job 6661462 at 5eaf0c5; 57 and 78 before the wider correction, job 6661283) (`slurm/h37-mutations.py`: each break replaces one piece of one source
  file and a test must fail), both in job 6661283 from the fresh clone. The first full run, at 5dab5dd (jobs 6661196 to 6661202 and 6661226 to 6661231), had 56
  tests and 74 breaks; the earlier break runs found six survivors, which were fixed by a test or by removing a guard that could change nothing.
- Two runs from fresh clones at f9f8fdd, jobs 6661281 to 6661287 and 6661288 to 6661293: every output file, the 22 tRNA alignments, the features, the
  contrasts and all tables, is identical byte for byte (`diff -r` of the two output folders). The first runs, at 5dab5dd, gave the same files, apart from
  the per-gene summary, which that commit did not write.
- The rest of the repository's tests, in the same job (607 passed, 183 skipped for want of inputs or tools, 2 failed): `tests/test_audit.py` fails
  on its comparison of the records' annotation with the build's changelog (the African grey parrot's patched record) and fails the same way at the
  parent commit 0d65412 on the same build (job 6661243), so it is this build and not this change; `tests/test_cybrid_swap.py` needs `torch`,
  which the environment of this job does not install. No existing module changed except `cli.py` (new commands) and `paths.py` (one folder); no earlier
  table is written by the new code, and no existing break script targets those two files.

## What is not assessed

- Part B and part C. Part C's precondition, that the rate pipeline reproduces H7 on the protein genes, is met; nothing of part C was started.
- Records with no annotation of a gene, the 94 Leu and 78 Ser copies with no isoacceptor, and the 20 odd-length features, as under "What each species carries".
- The reptiles' rate against a protein median (no median exists, as above); no cell was affected.
- The `all` group of H7's tables is not used by the bar; its RNA rows are in `rna_rates.tsv`.

## Files

Committed here: `gene_census.tsv`, `alignments.tsv`, `rrna_check.tsv`, `unresolved.tsv`; `rna_tests.<Group>.tsv` (every regression: gene, feature, model, set,
n, informative columns, lambda, slope, se, t, df, p, effect r, Holm p, q); `rna_rates.tsv` and `protein_rates.tsv` (H7's table format for the RNA genes and the 13
protein genes, on both sets, with the Holm p and the effect r); `rna_verdict.tsv` (one row per cell: the entry's own verdict
and, in the columns ending `_all`, the correction over all 768 cells and `verdict_all`); `rna_gene_summary.tsv`; `h7_control.<build>.json`; `sha256sums.txt`
(every output of the run, the 22 tRNA alignments, `rna_features.tsv` and the contrast tables included, which are not committed).

## Reproduce

From a fresh clone of the commit named in the entry, with the build, the tree, RefSeq release 237 and the CONTIG records in `data/input`, MAFFT v7.526 on PATH:
`slurm/h37-prepare.sbatch`, then `h37-genes`, `h37-control` (array of two), `h37-rates` and `h37-tests` (array of four), `h37-verdict`; `h37-check` runs the tests and the
breaks. The commands are in the repository's README.
