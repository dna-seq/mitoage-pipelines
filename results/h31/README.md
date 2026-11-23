# H31: the cybrids and the human variants by KronGenome's masked likelihoods

From the cluster. The likelihood tables are Slurm job 6591903 (one RTX 2080 Ti of an idle node, 10 minutes with the
second run; 12 s and 1.35 GiB a table), run from commit 0e5f0a4 with `slurm/h31-kron-table.sbatch`; the two runs gave
the same bytes for all six tables. The readouts (`mitoport kron-cybrids`, `mitoport kron-variants`) ran on the login
node from a fresh clone at commit 637762f, where the tests (47 passed with those of the nucleotide baselines and H32)
and the breaks of each guard (19 of 19 caught) ran too. An earlier job on a CPU node (6591902) failed as it had to:
the model's core runs on a GPU only. The rules were written before the code: `docs/HYPOTHESES.md`, H31.

| file | one row is |
| --- | --- |
| `tables/<accession>.mito.tsv` | one place of a host's mitochondrial genome: the log-probability of each base there when it is hidden (context: genetic code 2) |
| `kron_scores.tsv` | one (host, donor, complex): the sum and the mean of the drops over the complex's genes, how many places were scored and how many columns had a gap, and the pair's rRNA, tRNA and control-region sums |
| `h31_calibration.tsv`, `h31_summary.tsv`, `h31_verdict.tsv` | H16's readouts of those scores beside the baselines, and the rule |
| `kron_variants.tsv` | one human single-base variant: harmful or harmless, its kind of place, its effect in a protein gene, the drop, and how many bases the build's mammals hold there |
| `h31_variants_summary.tsv` | one area or comparison of the variant test, and the rule |
| `nuclear-arm/` | the same summaries with the nuclear context in place of genetic code 2 |

SHA-256 of the tables not kept here: NC_005089.1.nuclear 6922b04b…, NC_005089.1.shuffled 06d5c8ee…,
NC_012920.1.nuclear 06c795dd…, NC_012920.1.shuffled b27d192c… (the job's log has them whole).

Result: neither rule is met. Cybrids: `kron_cds` T1 −0.01, T2 0.59, T3 3 of 6; `kron_cds_mean` −0.05, 0.42, 2;
the bars were −0.63, 0.83 and 5. Variants: the area over all 95 harmful and 1,576 harmless variants is 0.88; on the
1,251 the conservation baseline covers it is 0.95 against 0.91, and the 95% interval of the difference, −0.015 to
0.095, holds 0. Reported beside: amino-acid-changing variants 0.93 against 0.87; tRNA genes 0.77, with no baseline.

## The shuffled genome, 2026-11-23

`arms-20261123/shuffled/` and `arms-20261123/mito/`: the readouts of Slurm job 6596363 (CPU, one minute), from a
fresh clone at commit c8b5b0c, on the tables of job 6591903. The shuffled genome's table is read at the record's own
places with the record's own bases (`--arm shuffled`); `mito` is the arm with genetic code 2 run beside it, on the 21
stated activities and the corrected baselines. Tests 21 passed, breaks 25 of 25 caught. Variants: area 0.364 on the
shuffled genome against 0.883 on the real one. Cybrids: `kron_cds` T1 0.06, T2 0.42 against −0.04, 0.59. The rule
was written before the scores: `docs/HYPOTHESES.md`, H31.
