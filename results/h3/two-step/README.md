# H3's two-step check of Coevol

From haumea, job 6490192, run from commit d6581f1 with `slurm/branch-rates.sbatch` on 2026-11-15. The inputs
are `coevol-prep`'s, byte-identical to the published Coevol run's (mammals.ali, .tree and .lht), from
MitoAge build `longevity-genie/mitoage` at c525adb and the dated tree from `dna-seq/mitoage` 6efdae4.

- `branches.tsv`: `longevity branch-rates`. One row per branch of HyPhy's MG94 fit with α and β per branch,
  with dS, dN, dN/dS, duration (Myr) and dS per Myr.
- `model.json`: S and N of that fit, its log likelihood and parameter count.
- `pgls.tsv`: `longevity branch-pgls`, the terminal branches against lifespan and mass.
- `rerun/`: the same two steps run again in the same job. They are not byte-identical (HyPhy's threaded
  optimiser), and agree to three significant figures.
- `slurm-6490192.log.txt`: the job's log.

The reading was fixed in docs/HYPOTHESES.md (H3) before the run. Its positive control (dN/dS rising with
body mass) fails, so the lifespan rows are not read.
