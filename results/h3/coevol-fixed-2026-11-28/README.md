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
of one that needs more time.

## The mechanism, corrected the same day (2026-11-28)

The paragraph that stood here gave the mechanism as: holding the ages fixed removes the parameter that would
otherwise absorb variation in rates, so the chain has nowhere to put it. That was offered as plausible, and an
outside review of the setup — asked for by the maintainer — found something far more ordinary sitting in our own
files. It is checked here against the data rather than taken on the reviewer's word.

**The dated tree this form was given has nine internal branches of exactly zero duration.** Verified two ways: by
counting the lengths in the Newick file (9 of its 104 non-root internal branches are `0`, and no tip is), and by
Coevol itself, which printed `warning : non strictly positive branch length : 0` into our own chain logs
**3,672 times in chain 1 and 3,744 times in chain 2**, from the first leg onward. Nobody read those lines.

With the ages held, a branch of zero duration still carries substitutions, so the rate on it is being asked to be
infinite. That is a degenerate parameterisation, and it is exactly the kind of thing that stops a chain mixing on
anything involving tree length while leaving other parameters looking healthy — which is the pattern we have.

**The ages-estimated form is not affected, and its own warnings are a red herring**: its tree carries no dates at
all, so all 104 internal branches are zero and Coevol warns about every one of them (82,317 times). There the ages
are sampled, zero is only a starting value, and that run converged. Zero durations are harmless when the ages are
free and degenerate when they are held — which is why this showed up in one half and not the other.

What remains unverified: the review also reported that a single terminal branch, *Erinaceus concolor*, accounts for
most of the late rise in the length trace. That branch is in the tree at 65.17 time units, long but not extreme
(the longest tip is 83.51 and the median 35.39), and the claim rests on a reconstruction of per-branch
contributions that has not been reproduced here. It is recorded as the reviewer's, not as ours.

**What this changes.** The chains stay unreadable and the stop stands: nothing can be read out of what was run.
But the reason is no longer "fixing the ages is inherently hard here". It is a defect in the input that can be
repaired, so declaring this half permanently unanswerable would be wrong, and the write-up must not say it. What
the analysis needs is a dated tree with strictly positive internal durations, and then a rerun.

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

It also does not say the question cannot be answered. With the mechanism corrected above, the honest statement is
narrower and more useful: *this run* was given a defective time tree and cannot be read, so the analysis has not
yet been done properly rather than having been done and failed.

## Where the review came from, and what it was

The maintainer asked on 2026-11-28 for an outside review of the Coevol setup. He asked for Astra-6.2; **that model
is refused for this account** ("not supported when using Codex with a ChatGPT account"), so the review ran on the
Codex model the account does serve, GPT-6, and the substitution was reported to him rather than passed over. It was
given the job script, the preparation code, the driver and both runs' convergence output, and no data files; it
read further files of its own accord, read-only, within this project and the Coevol source. Its claims are treated
as leads to check, not as results: the zero-duration branches are ours now because they were checked against the
tree and the logs, and the hedgehog branch is not.
