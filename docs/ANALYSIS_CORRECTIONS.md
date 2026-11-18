# Analysis corrections, 2026-11-18

The code review found five cases in which missing data or ties changed the reported answer. Historical result
tables and logs remain records of their original runs; the affected outputs need to be regenerated with these fixes.

- **H7 sign tests:** pairs with a zero rate difference now leave the sign-test denominator, as pairs with a zero
  lifespan difference already did. They stay in the regression. An entirely tied rate has no sign-test p-value.
- **H7 and adaptation outgroup contrasts (H21–H26):** the two sisters' distances to their outgroup now use one set of
  codon columns, resolved and non-stop in all three species. Independent pairwise deletion could invent a rate
  contrast when the sisters differed only in missing sequence. The existing minimum of 30 codons applies to this
  shared set. Pairwise distances for tests without an outgroup are unchanged.
- **C2 coevolution readouts:** a missing or infinite observed contrast is not assessed. Only finite null contrasts
  enter the empirical test, which still requires at least 100 valid draws; invalid draws are counted. Comparisons
  with NaN previously could produce the smallest possible empirical p-value.
- **Cybrid T4 readouts (H16–H18):** a total requires a finite score for every required complex, and every respiring
  comparator must have a total. H16/H17 require complexes I, III, IV and V. H18 explicitly requires its available
  host structures: four human complexes and three mouse complexes (I, III, IV); its baselines use the same sets.
- **H12 structural regressions:** an entirely unassessed run writes the usual result columns with missing estimates
  and p-values, instead of failing when it sorts a nonexistent p-value column. Printed unassessed rows are also
  selected after sorting.

`tests/test_pipeline_corrections.py` exercises these cases without downloads, GPUs or external executables. The
existing adaptation test now checks the two shared-mask distances against their output rows.

Regenerate `longevity rates`, `adapt pairs`, and `adapt longevity-pairs` for the changed rate calculations. Regenerate
`adapt coupling` and `adapt longevity-coupling` from their existing branch lengths, then their readouts, for C2.
Regenerate the cybrid readouts from their existing score tables for T4; the underlying coevolution, ESM and energy
calculations do not need to be repeated for this correction. Existing RELAX and aBSREL fits are unaffected.

The [corrected tables and validation record](../results/review-20261118/README.md) contain the remote rerun on build
`0f4e970`: H7, H12, H26 and the available cybrid readouts. The full test suite there passed 331 tests and skipped
142 tests whose optional inputs or dependencies were absent; the log lists the reasons. H21-H25 production C2 inputs
were absent, and H18 energy tasks were still incomplete, so those limits are recorded with the results.
