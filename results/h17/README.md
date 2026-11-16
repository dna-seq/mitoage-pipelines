# H17: the published cybrids scored by co-evolved contact pairs and by ESM-2

From haumea, Slurm jobs 6493605 (CPU: `mitoport cybrid`, `cybrid-coev`, the tests, 11 of 11 breaks caught, the rerun
identical) and 6493606 (GPU: `mitoport cybrid-esm`, `cybrid-readout`), both run from commit 329d9a6 on build 0f4e970,
with the contact maps and vertebrate partner alignments of H14's job. The logs are beside the tables. What was fixed in
advance and what was read: `docs/HYPOTHESES.md`, H17.

| file | one row is |
| --- | --- |
| `coev_scores.tsv` | one (host, donor, complex): the loss over contact pairs, over control pairs, their difference, the sensitivities, and the pairs scored, changed and unscored |
| `coev_pairs.tsv` | one contact pair with its control: the columns, the species and the effective sequences behind it |
| `esm_scores.tsv` | one (host, donor, complex): the ESM-2 sums alone and with the partner |
| `esm_llr.tsv` | one differing residue in one context (alone, or beside a nuclear subunit): the model's log-likelihood ratio |
| `h17_calibration.tsv`, `h17_summary.tsv` | H16's calibration and readouts T1 to T4 for these scores and H16's baselines |

`esm_checkpoint.txt` names the model; its commit was not reported by the library (the weights are those H12's job cached).
