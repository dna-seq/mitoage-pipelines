# H36, step 1: which nuclear subunits face the donors' mitochondrial differences

Slurm job 6616526 (CPU, one and a half minutes), run from a fresh clone at commit 3fb9e76 (`mitoport where`, then
`mitoport nuclear-where`, run twice: the three `nuclear_*` tables have the same bytes both times). Tests 22 passed and
the breaks 21 of 21 caught, on the same commit. The rules were written before the run: `docs/HYPOTHESES.md`, H36, with
the changes of 2026-11-24.

| file | one row is |
| --- | --- |
| `nuclear_subunits.tsv` | one nuclear subunit of one complex: its area for telling reduced from not reduced cybrids, the range over donors left out, the mean of the area under the permutation null, the subunit's and the family-wise p, whether it is named, whether that changes if the pairs the papers disagree on are coded reduced, its profile group |
| `nuclear_exposure.tsv` | one (host, donor, complex): how many residues of each nuclear subunit face a donor difference |
| `nuclear_residues.tsv` | one nuclear residue: in how many reduced and not reduced pairs it faces a difference (descriptive only) |
| `where_variants.tsv` | harmful and harmless human variants by what their site is (conserved, buried, at a contact) |
| `where_cybrid_sites.tsv`, `where_cybrid_summary.tsv` | one residue where a donor differs from its host with the description of the site; the summary per (host, donor, complex) |

Result: no subunit is named. See the entry.
