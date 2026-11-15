# H11 results from haumea

`mitoage longevity generation-time`, both maturity modes, run as Slurm job 6487335 on haumea
(`slurm/generation-time.sbatch`) from commit 858974b of this repository, in the clone
`/data/uq2347/runs/mitoage/pipelines-main`, on longevity-genie/mitoage at c525adb (species.parquet
794c421dfcdc) and the dated tree of `longevity dated-tree` (dated.nwk 5fc36105b0d1384d).

- `generation_time.{female,any}.tsv`: every fit (class, trait set, tree, model, term, feature, slope, p, q).
- `generation_time.{female,any}.classes.tsv`: per class and trait set, the species, each trait's correlation
  with lifespan, the variance inflation, and the loadings of the two slow-fast axes.
- `slurm-6487335.log.txt`: the job's log, with the commit, the tests (14 passed) and the rerun check
  (identical bytes for both modes).

What they show is in `docs/HYPOTHESES.md`, H11.
