# The manuscript's analyses on the build with the parrot's ND6

The same job as `results/euring/` (`slurm/manuscript-rerun.sbatch`, code at commit 2ee27a1 of this repository, so the
code is unchanged and every difference comes from the build), run once more on haumea on 2026-11-15:

| Run | Slurm job | Build | species.parquet | dated tree |
|---|---|---|---|---|
| after (`results/euring/`) | 6489877 | dna-seq/mitoage fb54cec = longevity-genie/mitoage f75bdc9 | 4e4aab66ae53 | fb54cec |
| nd6 (this folder) | 6490079 | dna-seq/mitoage 03c0949 = longevity-genie/mitoage 0f4e970 | 267a16311cca | a9819e2 (byte-identical to fb54cec's) |

The build is the EURING build plus the site repository's gene-annotation patch (its `docs/CHANGELOG.md`, 2026-11-15):
the African grey parrot (*Psittacus erithacus*, NC_027847.1) gains the ND6 its record leaves unannotated, and the ND6
alignment is redone with it (2,324 sequences). The data were fetched with
`mitoage download build --revision 0f4e970554b2fe0265e718321f422347293dde17`; the log records the digests, which equal
`public/data/` at 03c0949. The job ran in a detached worktree of 2ee27a1, so the log's first line names no commit.

- `*.tsv`: the tables the manuscript reads.
- `comparison.txt`: `slurm/compare-runs.py <after> <nd6>`, table by table and class by class.
- `slurm-6490079.log.txt`: the job's log.

**What the parrot changes.** With all 13 proteins it enters every bird set: birds 381 to 382.

- Unchanged, byte for byte: `aa_pgls` in mammals, fishes and reptiles; `fourfold.tsv`; every `generation_time`
  table (the parrot has no maturity age in the sets); the mammal, fish and reptile rows of `motif_pgls` and
  `structure_pgls`, so the reptile contact result (q = 0.046) is the one `results/euring/` gives.
- H2: mammals 5,262 tests, the same two columns pass with the same q (ND5 200, 0.0050; COX1 419, 0.049); 344 rows
  differ only in the ND6 columns, which the realignment renumbered. Birds 2,742 tests (was 2,738), none passes
  (smallest p 0.00066).
- H1, birds (382): nothing at q < 0.05 (smallest p 0.018).
- H7: 1,372 species with a tip, 340 cherries (102 bird pairs; was 339 and 101), dS > 1 in 24.1% of the concatenated
  pairs (was 24.2%), 655 tests. The same four pass: all species ND6 −0.197 (q = 0.0072), birds CYTB per million years
  −0.022 (q = 0.0085), mammals ND5 and ND2 unchanged. Birds, concatenated, per million years: p = 0.071, q = 0.27
  (96 pairs).
- H8: 132 tests, none passes; birds' smallest q 0.39 (was 0.41).
- H12, birds (382): contact minus control −0.016 (95% −0.039 to +0.007, q = 0.17; was −0.015, q = 0.20). Hydropathy:
  none under the primary model (COX2 q = 0.085); COX2 passes under the gene's own composition (q = 0.045, was 0.056)
  and with mass only (q = 0.030), as it did before the EURING birds (q = 0.041 and 0.027).
