# H16: the published cybrids against mismatch at the mito-nuclear contacts

Written by `uv run mitoage mitoport cybrid` on a laptop on 2026-11-15 (about 10 seconds), from the code of the commit
that adds this folder, and copied here from `data/output/mitoport/cybrid/`. What the tables mean and what was read
from them: `docs/HYPOTHESES.md`, H16.

| file | one row is |
| --- | --- |
| `positions.tsv` | one human position of a complex's mtDNA-encoded genes: resolved, touching a nuclear subunit (with how many nuclear residues and which chains), touching another mtDNA-encoded chain |
| `sequences.tsv` | one (species, gene): the accession its protein was read from, and whether the record is partial |
| `scores.tsv` | one (host, donor, complex): differing positions at contacts, at the other resolved positions and over all positions, with the counts compared |
| `calibration.tsv` | one measured (host, donor, complex): reduced, not reduced or the sources disagree, the stated activity where there is one, and the scores |
| `summary.tsv` | one readout of T1 to T4 (and the exploratory T2b) |

Inputs (SHA-256): `data/input/cybrid/records.gb` 0f7ffe14c23ab28b6f0a99f7961930de5b6e25659415aaecfb41883b26fbde6a;
`5XTD.cif` 7bee4569638b05a2b23709f17b51b28ed56c50b0d2df19d86028a490ad29ac32; `5XTE.cif` 2daafe97179293a3de6e75d562178b905f8f1c053b8e39e2b9b15c46b28bc6cb;
`5Z62.cif` c9401e6cbf48cb0fd0cc8df22f233fd2d894c295be9e144b03781324cb02adfa; `8H9S.cif` 8125baa6e3319d70ff30f7685b0be9690364a06b8868c9f29c5892fb30c45f10.
The measured values are `src/mitoage_pipelines/mitoport/cybrid_measurements.tsv`.
