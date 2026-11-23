# H34, part A: the predictors that already exist, on H29's variants

Run on the login node from a fresh clone at commit 832eca5 (a table join, seconds): `download mitimpact`
(MitImpact 4, archive and table pinned by SHA-256), then `mitoport variant-compare` on H29's corrected
`variant_scores.tsv` and H31's corrected `kron_variants.tsv` (`results/review-20261122/`). Tests 7 passed; breaks
10 of 10 caught. Reported, claimed as nothing: `docs/HYPOTHESES.md`, H34.

`h34_existing.tsv`: one row per score, with the variants it covers, its area under the ROC curve for harmful
against harmless, the area of the site's variability among the build's mammals on the same variants, the
difference with its 95% bootstrap interval, and what the score knew of labels like ours.

| score | harmful, harmless | area | baseline there | difference (interval) | labels |
| --- | --- | --- | --- | --- | --- |
| APOGEE 2 | 39, 328 | 1.000 | 0.919 | 0.081 (0.030 to 0.146) | trained on such variants of these genes: not a test |
| AlphaMissense | 39, 325 | 0.979 | 0.918 | 0.061 (0.015 to 0.117) | calibrated with population frequency |
| ESM1b | 39, 313 | 0.937 | 0.919 | 0.018 (-0.057 to 0.095) | none |
| KronGenome (H31) | 39, 328 | 0.923 | 0.919 | 0.004 (-0.058 to 0.072) | none |
| PhyloP, 100 vertebrates | 39, 328 | 0.838 | 0.919 | -0.081 (-0.155 to -0.010) | none |
| gnomAD local constraint | 39, 328 | 0.735 | 0.919 | -0.184 (-0.306 to -0.071) | built from population variation |
| StaB-ddG (H29) | 38, 314 | 0.629 | 0.915 | -0.286 (-0.422 to -0.150) | none |

Reading: the two scores that saw no labels and read sequence alone (ESM1b, KronGenome) are level with counting
the amino acids a site holds among mammals. Only the two that were fitted to labels like ours are above it, and
that is what such fitting gives. The bar a new model has to clear is therefore the baseline's 0.92, as H34 says.

# H34, part B: ProteinMPNN's membrane checkpoint on H29's variants

From the cluster's CPU nodes, code at commit b9e62ed (`mitoport variant-mpnn`, four complexes as four tasks, eight
cores each): first run, Slurm jobs 6596441 and 6598251; second run into another folder, jobs 6596442 and 6598252.
Four of the eight tasks of the first submission died at start on a missing torch file, because each task's `uv run`
re-synced the environment the others were importing from; they were run again with `UV_NO_SYNC=1` (6598251,
6598252) and the job file was fixed after (65406c6). All four tables of the two runs are identical byte for byte.
About 40 minutes of eight cores per run, most of it complex I. The collection, the readout, the tests, the breaks and
the check against the repository's own `score.py`: job 6600204 from commit 04b161d on the same tables (tests 47
passed, 1 skipped; breaks 22 of 22 caught). Two checks before it gave the same two tables and stopped on a break
that survived, each time a gap in a test and not in the code: 6598253 (variants that cannot be placed were tested
on complex IV, where every variant is placed; moved to complex V, 55bb14d) and 6600120 (a relative tolerance let
one seed's order pass for every seed; made absolute, 04b161d). The rule was written before the code: `docs/HYPOTHESES.md`,
H34.

| file | one row is |
| --- | --- |
| `mpnn_scores.tsv` | one variant in one gene: harmful or harmless, its place or why it has none, the score of each arm (the mean over seeds 0 to 2 and the copies of the chain), BLOSUM62, the site's variability among the build's mammals, contact or not |
| `h34_mpnn_summary.tsv` | one area or comparison, and the rule |

The score of each arm is the log-probability of the reference residue minus that of the new one, the position
decoded after every other residue of H29's system. `memb` is the global-label membrane checkpoint with the label
transmembrane (the score); `memb0` the same with the label soluble, `pmpnn` ProteinMPNN, `smpnn` the soluble-only
weights (the file StaB-ddG's weights started from).

Result: does not count. On the 38 harmful and 314 harmless variants every arm scores (the same variants as H29), the
area under the ROC curve is 0.583 for the membrane checkpoint, against 0.915 for the number of amino acids the site
holds among mammals: difference -0.332, 95% interval -0.472 to -0.188. Reported, deciding nothing: the label soluble
0.583, ProteinMPNN 0.581, the soluble weights 0.625, BLOSUM62 0.635; the membrane checkpoint minus the soluble
weights -0.041 (-0.095 to 0.003); at contacts with nuclear subunits 0.80 (8 harmful, 104 harmless; the baseline
0.89), elsewhere 0.54 (baseline 0.92); complex I 0.53 (29, 161; baseline 0.93), complex V 0.80 (9, 55; baseline
0.92); complexes III and IV hold no harmful variant.

Our reading against the repository's own `score.py` (job 6600204, `slurm/h34-mpnn-tool-check.py`): on its example
1BC8 and on COX2 with its smallest nuclear neighbour from 5Z62, each of the four arms, three draws, every place the
repository decoded last (279 and 810 of 279 and 810), the largest difference of a log-probability is 2.9e-6.

Reading: the membrane-trained checkpoint does no better than the soluble weights it was meant to improve on (and
the label that tells it the protein is a membrane protein changes almost nothing: Spearman 0.96 between the two
labels), and all four arms sit with BLOSUM62, far below counting the amino acids a site holds among mammals. The
contacts' 0.80 is eight harmful variants and below the baseline there too. Among the highest scores are common
CYTB alleles where the reference sequence carries the rarer residue; among the lowest, disease variants of ND5
and ND6 at sites that never vary.
