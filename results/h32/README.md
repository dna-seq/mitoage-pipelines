# H32: the conserved sequence blocks of the control region against random windows

From the cluster's login node (seconds), `mitoage mitoport cybrid-blocks` run from a fresh clone at commit 7407801
on the records of `download cybrid`; the tests (9 passed) and the breaks of each guard (11 of 11 caught) ran on the
same commit. The rule was written before the code: `docs/HYPOTHESES.md`, H32.

| file | one row is |
| --- | --- |
| `h32_blocks.tsv` | one (host, donor): the differences within the three blocks, in each block and in the whole control region; `places` are the blocks' places within the host's control region, counted from its first base |
| `h32_summary.tsv` | one check: the blocks' value and the share of 10,000 random window sets that reach it |

Result: does not count. On the mouse host, 16% of random window sets in the control region reach the blocks' area
of 0.98 and 59% their Spearman of −0.69; the rule needed under 5% for both. On the human host the bonobo differs
from human in none of the blocks' 58 places and the gorilla in 22, both with respiring cybrids, and one orangutan,
without, in 20.
