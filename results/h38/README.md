# H38: genome architecture against maximum lifespan

Run of 2026-11-27. The bars, the confound and the positive control are those of entry H38 in `docs/HYPOTHESES.md` as committed in fb04c2d, unchanged; the
choices the entry left open are written in the entry (commit c7990d8) before any run. Code commit a5c46ab. Two runs from fresh clones of it (Slurm jobs
6665458 and 6665460) gave every output file identical byte for byte, and the same bytes as the first two runs, at c7990d8 (jobs 6665266 and 6665268). Tests and
breaks: job 6665459.

## In one paragraph

**No architecture measure counts.** Nine measures in four taxon groups gave 36 cells; 33 are null and 3 are not assessed (too few species off the most common value
to regress). The positive control reproduced the audit's flags with no exception. Under the correction over all 384 p-values of the run (Holm, with
Benjamini-Hochberg q beside it, in the same table as the raw p) no regression that takes the tree into account has a Holm-adjusted p below 0.16 (lambda estimated: 1.0).
The only cells with a Holm p below 0.05 are 12 of the raw regressions, which ignore the tree and are only reported for comparison: the non-coding length and share
and the repeat share of ray-finned fishes. They are not a finding: every one is null once the tree is in the model.

## Inputs

The build of H27 and H37 (Hugging Face `longevity-genie/mitoage` at 446dad61883937075ec0c7e1f3f8135f2962a1cf; species table SHA-256
fa967e33724af62dab5859eab48393404627e54aee951bcd7fe69895b86a6db1), the dated tree of H37 (`phylogeny.json` SHA-256
bfcd086948b3e6d44b485201e1cd1af4c536360eebc935fd851d9d3c2ae19e76), RefSeq release 237 (md5 checked by the code) with the nine CONTIG-only records,
Tandem Repeats Finder 4.09.1 (bioconda build h031d066_4; binary SHA-256 dc8e610afe6674f770639fa0bcf026e09ba1b51588de68b8710c3d3466859fa8, pinned in the code, which refuses another;
parameters `2 7 7 80 10 50 500`). `architecture_inputs.tsv` has the checksums of the tables. TRF's raw output for every repeat is `trf_repeats.tsv`.

## The positive control, before any lifespan number was read

`arch_control.tsv`, `arch_control.json`, `arch_control_class.tsv`. 2325 live species, none without its record. **Reproduced, no exception in any of the 9 checks.** At the full tier: all 54 flagged
records other than a reversed deposit are above zero from their family, the 1 reversed deposit is at zero, all 470 records of the 81 uniform families are at zero, and where a family's
arrangement is its class's the class distance equals the family distance (484 records). At the genes tier: 15, 4, 1640 records in 173 families, 1646. The arrangements of the measure are the
audit's `order.tsv` for all 2325 records at both tiers. The audit run on this build flags 55 + 19 records; H13's published flags (build 0f4e970, 54 + 18) are all among them, one more at each tier.

**A departure from the entry's wording.** The entry says flagged records must be above zero and uniform families at zero, while measure (a) is the distance from the class's commonest arrangement. Read
against the class, that cannot hold: a family whose own arrangement is not its class's is uniform and above zero from the class (full tier: 29 records of uniform families, 41 records in families unlike their class;
genes tier: 10 and 13), and 4 flagged records at the full tier and 1 at the genes tier sit at the class's arrangement, so are at zero from the class. The control was therefore applied against the
family, with the same function and strings, and the class-level cases are listed in `arch_control_class.tsv` and counted in the json. This reading was written into the entry before the run.

## What each measure covers

2325 live species. Unknown at the full tier: 1 (*Crocidura russula*, no tRNA annotated); at the genes tier: 0. No control-region feature annotated: 350 species (their `cr_length`, `cr_copies`
and `repeat_cr` are unknown); more than one: 125. Distance above zero from the class's arrangement: 512 species at the full tier (mammals 111, birds 53, ray-finned fishes 170, reptiles 72) and 110
at the genes tier (2, 18, 28, 5). Species in the primary cells (H1's default dataset, confidence-filtered; the cell's size is per measure, so smaller where the control region is unknown) are in the table.

## Results

Primary set (H27's confidence-filtered records), lambda estimated and lambda 1 (the regression equal to phylogenetically independent contrasts), annotation flag not in the model; the
Holm p is over all 384 p-values of the run. Slope in standard deviations of the measure per tenfold lifespan, body mass in the model. Species in brackets: those flagged by the audit.

| measure | group | species (flagged) | slope, lambda 1 | p, lambda est. | p, lambda 1 | Holm over run, lambda est. | Holm over run, lambda 1 | high-subset species | verdict |
|---|---|---|---|---|---|---|---|---|---|
| cr_copies | Actinopterygii | | | | | | | | not assessed: only 3 species off the most common value |
| cr_copies | Aves | 323 (11) | -0.809 | 0.23 | 0.00124 | 1 | 0.448 | 290 | null |
| cr_copies | Mammalia | 528 (12) | -0.275 | 0.986 | 0.512 | 1 | 1 |  | null |
| cr_copies | Reptilia | 83 (4) | 0.371 | 0.381 | 0.252 | 1 | 1 | 74 | null |
| cr_length | Actinopterygii | 196 (3) | -0.12 | 0.814 | 0.559 | 1 | 1 | 156 | null |
| cr_length | Aves | 323 (11) | -0.973 | 0.431 | 0.0131 | 1 | 1 | 290 | null |
| cr_length | Mammalia | 528 (12) | 0.35 | 0.407 | 0.148 | 1 | 1 | 461 | null |
| cr_length | Reptilia | 83 (4) | 0.283 | 0.437 | 0.437 | 1 | 1 | 74 | null |
| gene_overlap_bp | Actinopterygii | 232 (5) | 0.0591 | 0.175 | 0.887 | 1 | 1 | 184 | null |
| gene_overlap_bp | Aves | 369 (12) | 0.0392 | 0.801 | 0.801 | 1 | 1 | 332 | null |
| gene_overlap_bp | Mammalia | 599 (14) | -0.234 | 0.972 | 0.477 | 1 | 1 | 527 | null |
| gene_overlap_bp | Reptilia | 95 (4) | -0.141 | 0.681 | 0.681 | 1 | 1 | 85 | null |
| non_coding_bp | Actinopterygii | 232 (5) | 0.212 | 0.249 | 0.317 | 1 | 1 | 184 | null |
| non_coding_bp | Aves | 369 (12) | -0.308 | 0.47 | 0.12 | 1 | 1 | 332 | null |
| non_coding_bp | Mammalia | 599 (14) | 0.467 | 0.225 | 0.0545 | 1 | 1 | 527 | null |
| non_coding_bp | Reptilia | 95 (4) | -0.00719 | 0.897 | 0.986 | 1 | 1 | 85 | null |
| non_coding_frac | Actinopterygii | 232 (5) | 0.204 | 0.256 | 0.333 | 1 | 1 | 184 | null |
| non_coding_frac | Aves | 369 (12) | -0.464 | 0.1 | 0.0091 | 1 | 1 | 332 | null |
| non_coding_frac | Mammalia | 599 (14) | 0.573 | 0.15 | 0.0167 | 1 | 1 | 527 | null |
| non_coding_frac | Reptilia | 95 (4) | -0.0376 | 0.854 | 0.929 | 1 | 1 | 85 | null |
| rearr_full | Actinopterygii | 232 (5) | 0.249 | 0.458 | 0.273 | 1 | 1 | 184 | null |
| rearr_full | Aves | 369 (12) | -0.913 | 0.445 | 0.00189 | 1 | 0.68 | 332 | null |
| rearr_full | Mammalia | 598 (13) | 0.547 | 0.358 | 0.0963 | 1 | 1 | 527 | null |
| rearr_full | Reptilia | 95 (4) | 0.353 | 0.268 | 0.268 | 1 | 1 | 85 | null |
| rearr_genes | Actinopterygii | 232 (5) | -0.0084 | 0.914 | 0.914 | 1 | 1 |  | null |
| rearr_genes | Aves | 369 (12) | -0.395 | 0.474 | 0.165 | 1 | 1 | 332 | null |
| rearr_genes | Mammalia | | | | | | | | not assessed: only 2 species off the most common value |
| rearr_genes | Reptilia | | | | | | | | not assessed: only 3 species off the most common value |
| repeat_cr | Actinopterygii | 196 (3) | -0.291 | 0.572 | 0.119 | 1 | 1 | 156 | null |
| repeat_cr | Aves | 323 (11) | -0.188 | 0.178 | 0.396 | 1 | 1 | 290 | null |
| repeat_cr | Mammalia | 528 (12) | 0.48 | 0.0641 | 0.0247 | 1 | 1 | 461 | null |
| repeat_cr | Reptilia | 83 (4) | 0.895 | 0.147 | 0.0474 | 1 | 1 | 74 | null |
| repeat_genome | Actinopterygii | 232 (5) | -0.0808 | 0.259 | 0.78 | 1 | 1 | 184 | null |
| repeat_genome | Aves | 369 (12) | -0.274 | 0.0382 | 0.076 | 1 | 1 | 332 | null |
| repeat_genome | Mammalia | 599 (14) | 0.378 | 0.0525 | 0.0486 | 1 | 1 | 527 | null |
| repeat_genome | Reptilia | 95 (4) | 1.22 | 0.0587 | 0.00812 | 1 | 1 | 85 | null |

Over the 384 cells with a p-value, 128 per model: raw p below 0.05 in 6 with lambda estimated, 36 with lambda 1 and 40 raw; Holm below 0.05 in 0, 0 and 12; q below 0.05 in 0, 5 and 22. The five with
lambda 1 are the number of control-region features of birds (high set, with and without the flag; confidence-filtered set) and the full-tier distance of birds, with q from 0.0095 to 0.028 and Holm 0.16 to 0.68;
each fails the first condition because lambda estimated is null (p 0.23 and 0.44). The covariate (the flag, in 3 to 14 species per cell) moves no cell that would otherwise have counted, since none did; its slopes are in
`arch_verdict.tsv` and its rows in `arch_cells.tsv`.

## Limits

Unknown: whether the lambda-1 signals in birds (the number of control-region features, the full-tier distance) are annotation: the audit's flag cannot catch a poorly annotated record in a family it cannot judge (below three
records). The distance is zero in most species (a measure with 2 or 3 species off zero is not assessed), so the regressions on it are crude; no power calculation was made. The test species are H1's: a species whose record lacks a protein gene is
not in them. The control-region count is a count of annotated features, and 350 records have none. A circular genome is read as a line by TRF.

## Reproduce

From a fresh clone of the code commit, with the build, the tree, release 237, the CONTIG records and the taxonomy dump in `data/input` and TRF on PATH: `slurm/h38-run.sbatch` (audit steps, measures, control, tests), `slurm/h38-check.sbatch` (tests and breaks).
