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
