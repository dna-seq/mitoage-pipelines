# H3, Coevol with the dated ages held fixed: the fifth leg is finished and the chains are still unreadable (2026-11-28)

This is the half of H3's Coevol run that answers the question actually asked — with the internal node ages held
at the dated tree's values, the synonymous rate is a rate per unit of real time rather than rate and time mixed
together. It does **not** pass the convergence bar the entry fixed in advance, after five legs, and it is stopped
here. **Nothing in it is read, and nothing is claimed from it, in either direction.**

The other half, with the ages sampled rather than held, is finished, converged and already read
(`results/h3/coevol-2026-11-27/`): with body mass held, nothing about lifespan passes the threshold. That remains
the only answer H3's Coevol run has.

## What was run

- Both chains reached their full **15,001 points** after five legs. Leg 5 was Slurm jobs 6666187 (chain 1,
  completed after 17 h 34 m) and 6666188 (chain 2, 21 h 28 m), from the code commit the earlier legs used.
- `tracecomp -x 1500 chain1 chain2` on copies of the two traces, 13,501 points each after a burn-in of 1,500.
  The copies were read in a folder of their own, as the ages-estimated readout was, so nothing writes into a
  directory a job owns. The repository ignores `*.log`, so that run's output is not committed here any more than
  the ages-estimated readout's was; every number it printed for a parameter that matters is in the table below,
  and `sha256sums.txt` names the bytes of the two traces and of the log, so the run can be identified again.

## Why it is stopped

The bar, from H3's entry and fixed before any of this: **effective sample sizes above 100 and small relative
differences between the chains.** Neither holds.

| parameter | effective sample size | relative difference |
| --- | --- | --- |
| tree length | **4** | 0.24 |
| sigma 2×3 | **71** | 0.20 |
| sigma 2×4 | **71** | 0.18 |
| sigma 1×4 | **92** | 0.17 |
| sigma 3×4 | **98** | 0.27 |
| sigma 3×3 | 113 | 0.25 |
| sigma 4×4 | 113 | 0.26 |
| sigma 1×3 | 138 | 0.21 |

Four of the covariance entries are below 100 and the tree length is at 4. The relative differences of the
covariances run from 0.17 to 0.27, which is not small: these are exactly the entries the correlations would be
read from, so the part that fails is the part the result would come out of. Two entries do pass on both counts
(sigma 1×2 at 391 and 0.007, sigma 1×1 at 523 and 0.036), and a result cannot be assembled from the two that
converged while seven did not.

**The tree length's history is the reason this is a stop and not a wait.** It was 31 after the third leg, 4 after
the fourth, and 4 again now after a fifth leg of about 3,500 more points per chain. An effective sample size that
falls, or sits still, while a chain gets longer is the signature of a chain fighting itself on that parameter, not
of one that needs more time. There is a plausible mechanism, stated when the question was filed rather than
invented now: holding the ages fixed removes the very parameter that would otherwise absorb variation in rates,
so the chain has nowhere to put it.

## The decision and whose it was

Desk question #96 put this to the maintainer on 2026-11-27 with three options — stop, buy three more legs
(about 430 core-days), or change the setup to a calibrated prior and start over — and recommended stopping.
Its stated default on silence was: *let the fifth leg finish, check the number once, and if it is still under 100
stop there and write it up as unreadable.* The fifth leg finished, the number was checked once, it is still under
100, and **the default is applied**. The question is withdrawn as answered by its own default, not decided by me
on the merits: if the maintainer wants the other two options the chains are still on disk and can be continued.

## What this does and does not say

It says the ages-held form of this analysis **could not be read**, and why. It does **not** say that there is no
relationship between the synonymous rate and lifespan in real time, and the write-up must not let the two be
confused: an unreadable chain is not a negative result. H3's claim stands on the ages-estimated half alone, with
that limit named.
