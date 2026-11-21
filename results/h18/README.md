# H18: the published cybrids scored by an interface energy on the chimeric complexes

From the cluster. The energies are Slurm array job 6500884 (44 tasks, OpenMM on CPU, 16 cores each, 30 minutes to 12 hours
a task), run from commit 313f3e3. The collection, the readout, the tests (38 passed) and the breaks of each guard (16
of 16 caught) are job 6575325, run from commit d2c3646: the first check job (6500885, from 313f3e3) stopped at one
failed test, which three small jobs (6572017, 6573071, 6574125) traced to single-precision rounding on an unrelaxed
structure, not to the energy function. What was fixed in advance, what was added after the array and what was read:
`docs/HYPOTHESES.md`, H18.

| file | one row is |
| --- | --- |
| `energy_scores.tsv` | one (host, donor, complex): the structure, the chains and residues of the system, the substitutions made and not modelled, the interaction energy of reference and chimera with the spread of three starts, `ddE`, `ddE_sd`, the substitutions at the interface and `ddE_per_change` |
| `h18_calibration.tsv` | one (host, donor, complex) with a published measurement: the measurement beside H16's baselines and the energy scores |
| `h18_summary.tsv` | the readouts T1 to T4 for the energy scores and H16's baselines (T2b, within one complex, was added after the first run and is exploratory) |

The ten tasks run a second time from a new build (the control of 2026-11-21) are not here yet.
