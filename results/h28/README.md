# H28: the cybrids by StaB-ddG on the hosts' complexes

From the cluster, commit 644323f. The model ran on RTX 2080 Ti cards (11 GB) of idle nodes, one task a card:
first run, Slurm jobs 6591264, 6591265 and 6591267 to 6591272 (9.6 card-hours); second run, into another folder,
jobs 6591783, 6591785, 6591792, 6591793, 6591800, 6591807, 6591808 and 6591810. All 44 tables of the two runs are
identical byte for byte. A task took 2.6 minutes and 0.56 GiB for complex V, 4.6 and 0.92 for IV, 7.6 and 1.55 for
III, 32 to 34 minutes and 3.53 GiB for I. The collection, the readout, the comparison, the tests (50 passed, 4
skipped) and the breaks of each guard (26 of 26 caught) are job 6592268, in the same clone. The rule was written
before the code: `docs/HYPOTHESES.md`, H28.

| file | one row is |
| --- | --- |
| `stab_scores.tsv` | one (host, donor, complex): the model's scores, mean over the seeds, with each control's |
| `stab_by_seed.tsv` | the same for each seed |
| `stab_genes.tsv` | one ported protein of a pair: its own score (shown as no verdict: H29 failed) |
| `h28_calibration.tsv`, `h28_summary.tsv`, `h28_verdict.tsv` | H16's readouts of the scores beside the baselines, and the rule |

Result: does not count. `stab_all`: T1 -0.10, T2 0.56, T3 3 of 6; `stab_contact`: -0.24, 0.65, 3 of 6; the count of
amino acid differences `d_all`, the baseline to beat: -0.36, 0.74, 4. The same verdict on each of the three seeds.
