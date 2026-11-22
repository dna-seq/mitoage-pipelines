# The cybrid comparisons after the review of 2026-11-22: H16's count corrected

Slurm job 6592805, from a fresh clone at commit 9d82153, with the inputs and the model's tables of H28's run
(nothing of the model was run again): `mitoport cybrid`, `cybrid-stab-collect`, `cybrid-stab-readout`, the tests
(57 passed, 4 skipped) and the breaks of the new guards (5 of 5 caught).

What changed: the donor is read through the host and not through its own alignment to human, and the end residues
of a partial record are unknown. Fourteen (host, donor, complex) counts fall, most for mouse hosts: M. caroli
complex V 23 to 15 and complex I 194 to 188, M. spretus I 71 to 65, M. terricolor I 140 to 130, M. pahari I 237 to
226, Otomys I 325 to 314; chimpanzee and bonobo complex I lose the one false substitution.

What it does to the tests: `d_all` T1 -0.355 as before, T2 0.742 (was 0.744), T3 4 of 6 as before; `d_contact` T2
0.806 (was 0.800). H28's verdict table is the same on every row but the baseline's 0.742.

| file | what |
| --- | --- |
| `scores.tsv`, `calibration.tsv`, `summary.tsv` | H16's tables with the corrected count |
| `h28_summary.tsv`, `h28_verdict.tsv` | H28's readout against the corrected baseline |

Not made again yet: H17's ESM and co-evolution sums (they place the donor on human columns the same old way and
need the model run again), H31's and the nucleotide baselines' readouts against the corrected count.

## The Chinese hamster's record corrected (job 6592840, commit fd27213)

From a second fresh clone, with KronGenome's likelihood tables of H31's run (the model was not run again):
`mitoport cybrid`, `cybrid-dna`, `cybrid-blocks`, `kron-cybrids`, `kron-variants`, the tests (72 passed) and the
breaks (19 of 19, 22 of 22, 5 of 5).

Mouse against the Chinese hamster: `nt_trna_d` 220 (was 280), `kron_trna` 20.7 (was 254.7); its control region is
now scored (`nt_dloop_d` 279, `nt_csb_d` 17, `kron_dloop` 62.5), and H32's second check holds for it as for the
golden hamster (above every donor that respired, as in 54% of random window sets). The donor has no measured
activity: T1 and T2 of the nucleotide scores are the same to the last digit, and H31's and H32's verdicts stand.
H31's variant summary is unchanged byte for byte.

With the corrected count, `d_all` within one complex separates reduced from not reduced with an area of 0.60 (I),
0.57 (III) and 0.72 (IV).

| file | what |
| --- | --- |
| `dna_scores.tsv`, `dna_summary.tsv`, `dna_calibration.tsv`, `dna_correlations.tsv` | the nucleotide baselines, remade |
| `kron_scores.tsv`, `h31_summary.tsv`, `h31_verdict.tsv`, `h31_calibration.tsv` | H31's cybrid readout, remade |
| `h32_summary.tsv`, `h32_blocks.tsv` | H32, remade |

## Changes of a start codon are no amino acid changes (job 6592885, commit 498780b)

From a fresh clone, on the tasks' tables of H29's run and KronGenome's likelihood tables (no model run again):
`variant-stab-collect`, `variant-stab-readout`, `kron-variants`, the tests (29 passed) and the breaks (24 of 24,
22 of 22).

H29: five variants in a first codon leave (two harmful, three harmless); 38 harmful and 314 harmless are scored
with both baselines. The model's area is 0.629 (was 0.637), BLOSUM62 0.635, the site's variability 0.915; the
verdict stands. H31: among amino acid changes KronGenome is 0.922 against the baseline's 0.865 on 38 and 321
variants (was 0.925 against 0.871 on 40 and 323), the interval of the difference -0.009 to 0.130; the rule over
all variants is as it was.

| file | what |
| --- | --- |
| `variant_scores.tsv`, `h29_summary.tsv` | H29, collected and read out again |
| `kron_variants.tsv`, `h31_variants_summary.tsv` | H31's variant test, remade |

## H28's positive control (job 6592766, commit 21a56ba): `h28_control.tsv`

`slurm/h28-stab-control.py` on four complexes of the model's benchmark (1AO7, 1A4Y, 2G2U, 4G0N; 150 mutations; 184
more left out because the Protein Data Bank's entry numbers the residue otherwise than the benchmark's cleaned
file, all of 1LFD among them). One row is a mutation: the measured value, the authors' published prediction and
ours. Ours against the authors': Pearson -0.993 over all, -0.967 to -0.997 by complex. The sign is opposite
because the benchmark file keeps "negative is destabilising" and the authors' run script, like the pipeline,
turns it once. Against the measured values ours has a Spearman of 0.375 in size and the authors' 0.380. The
log's "mean absolute difference" compares numbers of opposite sign and says nothing; it is not used.

## Twenty-one stated activities (job 6593046, commit da06759): `activities-21/`

Eight activities of McKenzie et al. 2003 (doi 10.1093/molbev/msg132), which the text does not state, are read
from the bars of its Fig. 4 against the control bar: M. caroli 89, 96, 108 (complexes I, III, IV), M. terricolor
83, 98, 94, M. pahari 87 and 105 (I, III). The same reading returns the values the text states (rat 46, 37, 78;
Otomys 72, 44; M. pahari IV 59) and the 2000 paper's for M. spretus (78, 114, 85) within 1 point. The paper
calls none of the eight a decrease, so T2 and T3 are as before; T1 has 21 values where it had 13. This is a
change of the measurement table after every result here was known.

T1 (Spearman of activity against the score; negative is the prediction), 21 values, with the 13-value figure of
the corrected tables in brackets where it was computed: `d_all` -0.520 (-0.355), `d_contact` -0.501 (-0.195),
`f_all` -0.475 (-0.594); `stab_all` -0.307 (-0.105), `stab_contact` -0.335 (-0.242); `kron_cds` -0.036 (-0.006);
`nt_cds_d` -0.370 (-0.036), `nt_nonsyn_d` -0.509 (-0.338), `nt_syn_d` -0.358 (-0.008), `nt_rrna_d` -0.454
(-0.652), `nt_trna_d` -0.501 (-0.643), `nt_dloop_d` -0.672 (-0.669), `nt_csb_d` -0.693 on 18 (-0.686 on 10).
H28's and H31's verdicts stand: neither model beats the count on any test. H32: the blocks' Spearman is reached
by 7.6% of random window sets of the control region (it was 58.8%), their area by 16.2% as before; the rule
asks both under 5%, so it still does not count.
