# Making H3's Coevol run converge sooner

H3 asks whether dS and dN/dS of the mitochondrial proteins follow lifespan and body mass along the dated tree
(Coevol 1.6, 106 mammals, one per family). The first two chains had not mixed after 1,100 points: the
effective sample size for dN/dS was 7. This note covers what was measured on haumea on 2026-11-15, the
alternatives we looked at, and what we chose.

## Where the time goes

- **Per point, serial.** A saved point costs about 90–96 s of wall time. The published run saved 1,100
  points in 29.3 h, and a benchmark on the test partition gave 90 s. The 76.6 s in `.chrono` counts only the
  MCMC moves. The rest is the likelihood over all sites, which Coevol recomputes on one core for the trace.
- **Per point, MPI build** (`Makefile.mpi`, gcc 15.2.0, OpenMPI 5.0.8). Same chain, same node type, a copy
  of the published chain 1:

  | Ranks | Wall time per point |
  |---|---|
  | 8 | 64 s |
  | 16 | 68 s |
  | 32 | 66 s |

  Rank 0 runs every parameter move on the pooled sufficient statistics, serially, so more ranks do not help.
  The MPI build read the serial chain's `.param` and continued it: the log likelihood picked up at
  −570,874 after −570,824. It never exits on its own, because rank 0 returns and the other ranks wait.
  `coevol-resume --mpi` ends it once its files are consistent and unchanged.
- **Mixing is the real cost.** At 900 post-burn-in points the autocorrelation of dN/dS was about 130
  points. Lartillot and Delsuc (2012, Evolution) reached ESS above 500 for every covariance entry with
  11,000 points. That implies an autocorrelation of 20 at most.
- **The ages.** `coevol-prep` gave Coevol the dated tree's topology without its branch lengths, so Coevol
  sampled about a hundred node ages with no calibration. Moving ages jointly with rates is the slowest part
  of its sampler. `-fixtimes` holds the ages at the dated tree's. That is also the question H3 asks, "along
  the dated tree".

## Alternatives

| Option | Same model? | Faster per point? | Mixes better? | Cost to adopt |
|---|---|---|---|---|
| Coevol 1.6, MPI build, 8 ranks | yes | 1.5× | no | done (`--mpi 8`) |
| Coevol with `-fixtimes` on the dated ages | yes, with ages given | some (no time moves) | expected, to be measured | done (`--fixed-ages`); new chains |
| BayesCode `nodeomega` (Latrille, Lartillot; bioconda, v1.3.4) | nearly: ω and mutation rate as a Brownian process with traits, its own codon model and priors | unknown, likely similar | same move types (one substitution mapping, then 30 parameter sweeps); no reason to expect better | medium: new inputs, results not comparable one to one |
| RevBayes | would have to be written (it has multivariate BM tutorials, not ω with traits) | likely slower | unknown | high, unvalidated |
| BEAST 1.10 (HMC relaxed random walks, factor analysis) | no: traits, not branch ω jointly with traits | n/a | n/a | wrong model |
| GPU, Rust or BEAGLE reimplementation of Coevol | none found | | | |
| Two-step: ML dS and dN/dS per branch, then PGLS | an approximation | hours instead of days | no MCMC | done (`branch-rates`, `branch-pgls`) |

## What we chose

1. Keep Coevol. Nothing faster implements the same model, and BayesCode's sampler is the same design.
2. Resume the two published chains on the MPI build with 8 ranks. Each chain runs in 3-day legs with a
   checkpoint every 6 h, towards 15,000 points: about 10 days at 64 s per point. This matches Lartillot and
   Delsuc's 11,000.
3. Run two new chains with the ages held (`-fixtimes`) beside them. If they mix faster, as expected, they
   give the converged answer first, and they answer the question as H3 states it. The resumed chains then
   show what sampling the ages changes.
4. Run the two-step check now (minutes on 16 cores). It is a check of direction and size; it ignores the
   uncertainty of each branch's rates. Its reading was fixed in H3 before it ran.
