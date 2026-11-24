# H36, step 2 for the mouse host: the subunit swap on the mouse's own solved complexes (StaB-ddG; rat, M. caroli, M. pahari)

Tasks from commit 7aec214: Slurm array 6631178 (all 171 tasks completed), second run 6631179 into another folder (171 of 171).
The two runs' 342 task tables are identical byte for byte, and so are the eight tables below. Collection, the byte comparison, the
tests (126 passed) and the breaks of each guard (96 of 96 of the swap's, 60 of 60 of the donor orthologs' and the swap's protein reading)
at commit 369700c: job 6638790. CPU nodes only, 2 cores and 16 GB per task (the first submission, jobs 6630435 to 6630437 with 8 cores,
was cancelled before any task started: the shared CPU cap kept it waiting; the job file says 8 cores and the submission overrode it).
The collection's code differs from the tasks' only by a refusal in `unassessed_donors` and a simplified expression in the reduced-pairs
table, neither of which touches what a task computes. The rule: `docs/HYPOTHESES.md`, H36, "Step 2 for the mouse host, made exact",
written before any score. The tasks took 94 wall-clock hours at 2 cores (about 63 core-hours per donor; the estimate of the first run, 8
per donor, was too low for a host that differs from its donors at many residues).

| file | one row is |
| --- | --- |
| `swap_scores.tsv` | one model reading: (donor, complex, nuclear gene, kind B0 / B1 / null, draw, seed), the three terms and the predicted change of binding energy `ddG` (positive: weaker) |
| `swap_cases.tsv` | one (donor, complex, nuclear gene): scored or why not assessed, residues placed and not placed, substitutions, contact differences, null candidates, whether the null was run, the mouse and donor RefSeq proteins |
| `swap_restoration.tsv` | one scored (donor, complex, gene): B0, B1, restoration, the null's mean and 95th percentile, share of draws at or above, rank among the complex's genes |
| `swap_candidates.tsv` | one (complex, gene): in how many donors it is among the largest and above its null, and the candidate rule |
| `swap_reduced.tsv` | one (complex, gene): the same rule on the pairs stated as reduced alone (rat I, rat IV, M. pahari IV), and the count in the not-reduced pairs |
| `swap_sums.tsv`, `swap_spearman.tsv` | the restoration summed over the genes scored for every donor, beside the stated activity, and the Spearman correlation |
| `swap_unassessed.tsv` | the donors with a stated outcome and no RefSeq source (*M. spretus*, *M. terricolor*, *Otomys irroratus*), for each complex, which have no task |

## What was and was not assessed

- 171 tasks (rat, *M. caroli*, *M. pahari*; complexes I, III, IV; every nuclear gene of the complex): 123 scored and 48 not assessed.
  The 48 are the genes with no chain in the mouse structure that touches an mtDNA-encoded chain and is the mouse's protein for that
  gene: I NDUFA2, NDUFA5, NDUFC1, NDUFS1, NDUFS4, NDUFS6, NDUFV1, NDUFV2, NDUFV3 (they touch none); III UQCR10, UQCRC2, UQCRH
  (none); IV COX6A1, COX7A2, COX8A (the structure's chains there are other isoforms, 90% identity not reached; not checked which) and
  NDUFA4 (no chain). Three donors each.
- Of the 123 scored: 30 have a donor protein that changes no placed residue (restoration exactly 0, null not run); 52 have a null of
  50 draws; 39 have fewer than 3 mammals in the null's window (null not assessed) and 2 have no alignment of the gene.
- *M. spretus*, *M. terricolor* and *Otomys irroratus*: not assessed, no ortholog source (9 rows, `swap_unassessed.tsv`). The Chinese
  and Syrian hamsters have only a whole-cybrid outcome (no respiring cybrid) into the mouse and are not cases; there is no hamster host.
- Not run: part (a) of the step (embeddings of a structure model; it needs a GPU), the controls with the other weights, a residue's own
  contribution.

## Result

No candidate, by the entry's rule and by its reading on the reduced pairs. Restorations (B0 minus B1, kcal/mol): 36 above 0, 57 below,
30 exactly 0; the seed-to-seed spread has a median of 0.04 and a 90th percentile of 0.13. Of the 52 with a null, 10 are above its 95th
percentile; 6 of those are also among the largest of the complex (top 3 for I's 28 genes: a quarter rounded up, 7; III 3; IV 3) and no
gene has that in more than one of its three donors (two are needed).

Top three genes by restoration per donor (null: its 95th percentile and the share of draws at or above the real value; "no null" says
how many mammals fell in the window):

| complex, donor | first | second | third |
| --- | --- | --- | --- |
| I, rat (reduced) | NDUFA1 0.44 (no null, 0) | NDUFS2 0.18 (p95 0.18, share 0.24) | NDUFB3 0.12 (p95 0.70, 0.29) |
| I, *M. pahari* | NDUFB4 0.66 (no null, 2) | NDUFA11 0.42 (no null, 0) | NDUFS2 0.14 (no null, 0) |
| I, *M. caroli* | NDUFA11 0.64 (no null, 0) | NDUFS5 0.13 (no null, 0) | NDUFB8 0.10 (no null, 0) |
| III, rat (mixed) | UQCRQ 0.42 (p95 0.64, 0.27) | UQCR11 0.12 (p95 0.45, 0.55) | CYC1 -0.09 |
| III, *M. pahari* | UQCRC1 0.59 (p95 0.49, 0.02) | UQCRQ 0.22 (p95 0.60, 0.59) | UQCR11 0.00 |
| III, *M. caroli* | UQCRC1 0.08 (p95 0.50, 0.49) | CYC1 0.00 | UQCRFS1 0.00 |
| IV, rat (reduced) | COX4I1 0.66 (p95 0.63, 0.02) | COX5B 0.12 (p95 -0.06, 0.02) | COX6C 0.03 (p95 -0.07, 0.02) |
| IV, *M. pahari* (reduced) | COX5A 0.00 | COX5B 0.00 | COX6B1 0.00 |
| IV, *M. caroli* | COX6C 0.11 (p95 0.15, 0.18) | COX5A 0.00 | COX6B1 0.00 |

- *The rule (more than half of the three donors):* no gene is both among the largest and above its null in more than one donor. The
  six that are in one donor: *M. caroli* I NDUFB1 (0.09, null p95 -0.15), *M. pahari* I NDUFB9 (0.03, -0.07), *M. pahari* III UQCRC1
  (0.59 against 0.49, a pair stated as not reduced), rat IV COX4I1 (0.66 against 0.63), COX5B (0.12 against -0.06) and COX6C (0.03
  against -0.07). In the last five the null's 95th percentile is itself negative or the restoration is within the seed noise.
- *The reading on the reduced pairs:* complex I has one reduced donor (rat): no gene is above its null among the largest (NDUFS2's 0.18
  equals its null's 95th percentile, NDUFA1's 0.44 has no null); complex IV has two (rat, *M. pahari*): COX4I1, COX5B and COX6C are
  both in the rat only, one of two, which is not more than half, and *M. pahari* IV, a reduced pair (59%), has no restoration above 0 in
  any of its 7 genes (four proteins identical to the mouse's, COX6C -0.06, COX7B -0.10, COX4I1 -0.21); complex III has no reduced pair.
  The not-reduced pairs hold three of the six (NDUFB1, NDUFB9, UQCRC1).
- *Against the activities:* the restoration summed over the genes scored for every donor of the complex, against the stated percent:
  complex I (n 3: *M. caroli* 89, *M. pahari* 87, rat 46) Spearman 1.00; IV (n 3) 0.50, p 0.67; III has two pairs with a percentage;
  pooled over complexes (n 8) 0.38, p 0.35. n is 3 or 2 in the complexes, so this claims nothing.
- *Reading:* as for the apes, the swap does not name a nuclear subunit for the mouse host. The two largest restorations with a null
  (rat IV COX4I1 0.66 and *M. pahari* III UQCRC1 0.59) exceed its 95th percentile by 0.03 and 0.10 kcal/mol, with no draw of 50 as
  large (the least a 50-draw null can show), each in one donor of three, and the other genes above their null do so within the
  seed noise or over a null whose own percentile is negative; the complex I values of 0.4 to 0.7 have no null to be compared with. The
  reduced *M. pahari* IV pair is the clearest absence: no nuclear subunit of the 7 with a chain, swapped for its donor's, restores the
  fit the model predicts. That says nothing of the four genes with no chain in the structure (COX6A1, COX7A2, COX8A, NDUFA4), nor of a
  cause that is not the binding of a subunit to the mtDNA-encoded chains.
