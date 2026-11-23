# H21 to H25: adaptations (diadromy, the deep sea, marine and hibernating mammals, bats)

From the cluster, a fresh clone at commit ff05158. The first job (6587083: tests, breaks 25 of 25, traits, pairs,
fits, selection inputs, all twice with the same bytes), the selection array (6587084: RELAX and aBSREL, 92 tasks
completed, none failed; the mammal-wide hibernation tasks were left out, their table is not obtained), the coupling
job (6587349, from commit 59f739a) and the readout (6599963, `mitoage adapt readout`). The rules were written before
the run: `docs/HYPOTHESES.md`, H21 to H25.

| file | one row is |
| --- | --- |
| `readout.tsv` | one readout of one hypothesis: its status, estimate, p, Holm p within the hypothesis, and the reading the entry fixed |
| `pair_tests.tsv`, `pairs.tsv`, `bat_pairs.tsv` | the tests over one pair per origin, and the pairs |
| `l2.tsv`, `lifespan.tsv` | the fits of H24 and H25 |

Result: no hypothesis is supported. See each entry.
