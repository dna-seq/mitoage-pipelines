# H3, Coevol: the first readout of the converged ages-estimated chains (2026-11-27)

The two resumed Coevol chains of H3 (106 mammals, one per family; lifespan and body mass as continuous
traits; the ages sampled, not held) reached 15,002 points each after five legs and pass the convergence bar
the entry fixed in advance: effective sample sizes above 100 and small relative differences between the
chains. They are read here for the first time. The chains with the dated ages held (`coevol-fixed`) do **not**
pass it and are not read; see H3's entry.

- Code commit the chains ran from: `d6581f1` (`slurm/coevol-chain.sbatch`, Coevol's MPI build on 8 ranks).
- Legs: jobs 6490194–6490201 (legs 1–4 of the two chains) and 6666186 (leg 5, chain 1).
- Readout: job 6666586, `readcoevol -x 1500 1 15000 chain<n>` on copies of the chains, 13,500 points read
  after a burn-in of 1,500. The copies were read in a folder of their own so that nothing wrote into a
  directory a running job owns.

## Files

- `chain1.cov`, `chain2.cov`: `readcoevol`'s covariance matrices over dS, dN/dS, lifespan and mass, with the
  correlation coefficients, the partial correlations given the other three variables, and for each the
  posterior probability that the entry is positive.
- `tracecomp.log`: `tracecomp -x 1500 chain1 chain2` on the same copies. Every effective sample size is 165
  or more (tree length 165, dN/dS 339, the covariances 1,350 to 10,028) and every relative difference is 0.16
  or less.

## The numbers (chain 1 / chain 2, posterior probability of a positive value in brackets)

Correlations:

| pair | r | pp |
| --- | --- | --- |
| dS × lifespan | −0.493 / −0.496 | 0.0003 / 0.00007 |
| dN/dS × lifespan | +0.439 / +0.439 | 0.99 / 0.99 |
| dS × mass | −0.381 / −0.384 | 0.0039 / 0.0027 |
| dN/dS × mass | +0.212 / +0.212 | 0.86 / 0.86 |
| lifespan × mass | +0.715 / +0.715 | 1 / 1 |
| dS × dN/dS | −0.633 / −0.634 | 0.00015 / 0.00007 |

Partial correlations, given the other three variables:

| pair | r | pp |
| --- | --- | --- |
| dS × lifespan | −0.116 / −0.120 | 0.15 / 0.14 |
| dN/dS × lifespan | +0.287 / +0.284 | 0.92 / 0.92 |
| dS × mass | −0.162 / −0.163 | 0.063 / 0.063 |
| dN/dS × mass | −0.227 / −0.225 | 0.11 / 0.11 |
| lifespan × mass | +0.679 / +0.678 | 1 / 1 |
| dS × dN/dS | −0.536 / −0.537 | 0.00096 / 0.0013 |

## What it says

With body mass held, nothing about lifespan reaches the entry's bar of a posterior probability of 0.975. The
rise of dN/dS with lifespan is +0.29 at 0.92, and the raw fall of dS with lifespan (−0.49, pp 0.0003) goes
with body mass: as a partial correlation it is −0.12 at 0.15. These are the numbers the unconverged chains
gave at 1,100 and at 4,300 points, to the second decimal, now from chains that have mixed; more points will
not move them.

This is a result and no longer "not assessed", but only for the analysis in which Coevol samples the node
ages itself. H3's entry notes that with the ages held dS is a rate per unit of the dated tree's time, which is
the reading that answers the question as stated, and that those chains hinted at something else. They have not
mixed and are not read, so that reading stays open.
