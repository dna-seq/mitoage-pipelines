# H29: StaB-ddG on single human variants, harmful against harmless

From the cluster, commit 4b56a5f, on RTX 2080 Ti cards (11 GB) of an idle node, one variant a call: first run,
Slurm jobs 6591811 to 6591814 (7.1 card-hours); second run into another folder, jobs 6592223 to 6592226. All 19
tables of the two runs are identical byte for byte. The collection, the readout, the tests (31 passed) and the
breaks of each guard (19 of 19 caught) are job 6592227 (a first try, 6592222, stopped for a missing input table).
The rule was written before the code: `docs/HYPOTHESES.md`, H29.

| file | one row is |
| --- | --- |
| `variant_scores.tsv` | one variant in one gene: harmful or harmless, its place in the structure or why it has none, the model's scores with each control's, and the two baselines |
| `h29_summary.tsv` | one area or comparison, and the rule |

Result: does not count. On the 39 harmful and 316 harmless variants scored with both baselines the area under the
ROC curve is 0.64 for the model's stability term (`stab_fold`), 0.58 for its binding term, 0.62 for plain
ProteinMPNN and 0.63 for the stability-only weights; 0.64 for BLOSUM62 and 0.92 for the number of amino acids the
site holds among the build's mammals. The model minus that baseline: -0.28, 95% interval -0.42 to -0.14.
Reported and claimed as nothing: at contacts with nuclear subunits the model reaches 0.87 (9 harmful, 105
harmless), away from them 0.59.

Known errors, found by the review of 2026-11-22 and not corrected here: two changes of a start codon are scored as
ordinary amino acid changes (m.7587T>C, harmful; m.3308T>C, harmless). Too few to move an area.
