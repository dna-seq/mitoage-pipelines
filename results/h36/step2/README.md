# H36, step 2: the subunit swap (StaB-ddG, human host, the apes' own nuclear orthologs)

Tasks run from commit 2c3d50c: Slurm array 6618903 (46 of its first tasks died at start on a file error in the freshly
synced environment and were run again as array 6619266; all 216 tasks completed), second run 6618904 into another folder
(216 of 216). The two runs' 432 task tables are identical byte for byte, and so are the six tables below. Collection, the
byte comparison, the tests (60 passed) and the breaks of each guard (64 of 64 caught): job 6625589 at commit 03c5d36 (the
first collection, job 6619267 at 2c3d50c, read the rows of kind `null` as missing values and was discarded; the task tables
were not affected). CPU nodes only. The rule: `docs/HYPOTHESES.md`, H36, "Step 2 made exact", written before any score.

| file | one row is |
| --- | --- |
| `swap_scores.tsv` | one model reading: (donor, complex, nuclear gene, kind B0 / B1 / null, draw, seed), with the three terms and the predicted change of binding energy `ddG` (positive: weaker) |
| `swap_cases.tsv` | one (donor, complex, nuclear gene): scored or why not assessed, residues placed and not placed, substitutions, contact differences, null candidates, whether the null was run |
| `swap_restoration.tsv` | one scored (donor, complex, gene): B0, B1, restoration, the null's mean and 95th percentile, share of draws at or above, rank among the complex's genes |
| `swap_candidates.tsv` | one (complex, gene): in how many donors it is among the largest and above its null, and the candidate rule |
| `swap_sums.tsv`, `swap_spearman.tsv` | the restoration summed over the genes scored for every donor, beside the stated activity, and the Spearman correlation |

Result: no candidate. See the entry. Part (a) of the step (structure-model embeddings) was not run: it needs a GPU.
