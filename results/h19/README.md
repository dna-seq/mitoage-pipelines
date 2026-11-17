# H19: the coupling in birds and fishes, against distance from the mtDNA-encoded core, and its mitochondrial side

From haumea, Slurm job 6493878 (`slurm/h19-classes.sbatch`), run from commit c63661b on build 0f4e970 (species table
267a16311cca) with the dated tree of dna-seq/mitoage 03c0949: the 26 tests passed, 33 of 33 breaks were caught, the
50 alignments of H14's job came out identical, and a second run gave the same bytes for every table. The log is beside
the tables. What was fixed in advance and how it is read: `docs/HYPOTHESES.md`, H19.

| file | one row is |
| --- | --- |
| `class_readouts.tsv` | one readout (C1, G1-G3, M1-M2, O1-O3, S1-S2, T1-T3) in one class: the groups, their medians, the estimate, p and Holm p |
| `class_coupling.tsv` | one nuclear or control protein in one class against one set of mitochondrial genes: species, pairs, columns, divergence and every coupling statistic |
| `class_matrix.tsv` | one group of proteins in one class against one set of mitochondrial genes: the median coupling, descriptive |
| `core_distance.tsv` | one nuclear subunit: the median and minimum distance of its heavy atoms to the nearest mtDNA-encoded chain, in Å |
