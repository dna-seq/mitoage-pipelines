# H36, step 3: SaProt 650M on the joined mitochondrial and nuclear chains

The model reads one mitochondrial chain, a linker of 25 glycines and one nuclear chain of the host's own solved complex, with the structure
tokens of Foldseek, and gives the pseudo-log-likelihood (PLL, nats) of the mitochondrial residues that touch the nuclear chain. The rule, the
bars and the null were written before any score: `docs/HYPOTHESES.md`, H36, "Step 3 made exact". **Result: the model does not count** (two of
the four conditions fail; see the entry).

Tasks: code commit c0ccb79. Probe job 6655272 (one RTX 2080 Ti: SaProt in 32-bit floats, attention `sdpa`, 3.5 GiB peak, about 3 s a read of
30 places, a read repeated gives the same bytes). First run: human host array 6655308, mouse host array 6655425; second run: human host
6655426, mouse host 6655427; the shards that stopped at start (17 of 139 array elements, a file of the Python environment not found at start,
none during a task) were resubmitted as 6655589, 6655590, 6655591, 6655684, 6655685, 6655720 and 6655804 and skipped the tasks that had their
table. No job was requeued by the node's owner (0 restarts). Cards of the requeue partition only, never an A100. The two runs' task tables
(864 for the human host, 684 for the mouse host) and all collected tables are identical byte for byte, although the second run was on other nodes
than the first (the card model of those nodes was not recorded: unknown). Collection, comparison, tests and breaks: job 6655958 at commit
98a9d8f (two bugs of the collection, found when the human host was first collected, were fixed in 98a9d8f; the tasks' code is unchanged since
c0ccb79). Tests: 116 passed (this step's 39 and the swap's). Breaks: 83 of 83 caught (the first break run, on the tests as first written, let 13 survive; their tests were added in 7ce19d2 and 98a9d8f; one change, the nuclear chain's substitutions put on every copy, is equivalent on these structures, no gene having copies that differ in sequence, and is not run). The swap's own breaks were not run again.

| file | one row is |
| --- | --- |
| `human/`, `mouse/` `saprot_scores.tsv` | one read of a task: B0 (the host's own nuclear chain), B1 (the donor's ortholog), a null draw, a shuffled partner; `pll`, and `ddG` = minus the PLL (the column the step-2 rule reads) |
| `saprot_pairs.tsv` | the same per pair of chains (mitochondrial chain, nuclear chain), for B0, B1 and the shuffled partners |
| `saprot_cases.tsv` | one (donor, complex, gene): scored or why not assessed, placed residues, pairs, pairs too long, pairs without structure states |
| `saprot_restoration.tsv` | one scored task: restoration = PLL(B1) - PLL(B0), the null's mean and 95th percentile, rank in the complex |
| `saprot_candidates.tsv`, `saprot_reduced.tsv` | the step-2 candidate rule on all donors, and on the reduced pairs alone (mouse host) |
| `saprot_sums.tsv`, `saprot_spearman.tsv` | the restoration summed over the genes scored for every donor, beside the stated activity |
| `saprot_residues.tsv` | one changed residue of a nuclear chain read alone beside one mitochondrial chain: change of the PLL, distance, interface |
| `saprot_residue_contrast.tsv`, `saprot_residue_contrast_pooled.tsv` | the interface contrast per host, and for both hosts' rows together |
| `saprot_dependence.tsv` | PLL(B0) minus the mean PLL with the partner shuffled |
| `saprot_compare.tsv`, `saprot_compare_cases.tsv` | the bar: SaProt's restoration, StaB-ddG's (step 2) and step 1's count, over the same cases and genes |
| `mouse/saprot_unassessed.tsv` | donors with a stated outcome and no RefSeq source |

The pooled contrast is made by `cybrid_saprot.interface_contrast` over the two hosts' `saprot_residues.tsv` rows concatenated.

## What was and was not assessed

- Human host (apes), 216 tasks: 149 scored (as in step 2), 66 without an alignment of the gene, 1 without the donor's ortholog (the bonobo's
  NDUFA13). Of the 149: 88 have an ortholog that changes no placed residue (restoration exactly 0), 61 change a residue, 40 of those have a null of
  50 draws. 315 pairs read; 3 pairs (complex I) longer than 1,022 residues were not read.
- Mouse host (rat, *M. caroli*, *M. pahari*), 171 tasks: 120 scored, 48 with no chain of the mouse's own protein touching the mitochondrial chains
  (as in step 2) and 3 with no readable pair; of the 120, 92 change a residue and 51 have a null. 252 pairs read; 3 pairs too long. *M. spretus*,
  *M. terricolor* and *Otomys irroratus*: not assessed, no ortholog source (`mouse/saprot_unassessed.tsv`). No hamster host.
- Four chains of the mouse complex I (8OM1) have Foldseek tokens that do not match the chain as read; no pair was lost to it (`pairs_no_states`
  is 0 in every case).
- Not read: the model's attention.
