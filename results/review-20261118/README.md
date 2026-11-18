# Corrected analyses from the 2026-11-18 review

Run on the remote Linux host from the corrected source after baseline commit `1708f8e`, on MitoAge build `0f4e970`.
The species table SHA-256 is `267a16311cca579af2ca1cd004e99befe66e53e1ae7696d71ffb9f6b7be5dbed`.
`validation.json` records the exact source and codon-alignment hashes. [The corrections](../../docs/ANALYSIS_CORRECTIONS.md)
explain the changed calculations. Original run directories and historical tables were preserved.

- H7: all 655 result rows were regenerated. Combining shared-codon deletion with the corrected sign-test denominator
  changes 32 rows' nominal sign-test p<0.05 classification and two rows' nominal regression p<0.05 classification.
  No regression q<0.05 classification changes. Rows repeat classes, genes and normalisations; they are not independent
  discoveries. `sign_test_changes.tsv` gives the 32 changed sign-test rows.
- H12: all 192 structural-regression result rows were regenerated. Every p-value agrees with the previous run.
  The entirely unassessed case is covered separately by the regression test.
- H26: primary unsaturated dN/dS contrasts remain not assessed in mammals (6 pairs; minimum 8) and not detected in
  birds (41 pairs; Wilcoxon p=0.662777, previously 0.672121). The saturated sensitivity tests were also regenerated.
- H16/H17: readouts were regenerated from existing scores. All 28 H17 T4 values are unchanged on the complete data.
- H18: the snapshot contains only 4 of 44 energy tasks. The partial collection and readout completed, and all energy
  T4 totals remain not assessed. They are validation of missing-data handling, not a final H18 result; see `rerun.log`.
- H21-H25 adaptation rates and production C2 readouts were not rerun because their prepared traits/branch lengths
  are absent here. The new shared-mask adapter and finite-null guards have direct regression tests.

Validation on the remote checkout: 61 focused tests passed; the complete suite with H12's prepared data and the
external-tool environment on PATH finished with 331 passed and 142 skipped. `pytest.log` lists every skip reason.
The additional regression tests need no external programs or downloads.

`rerun.py` reads inputs through `MITOAGE_PIPELINES_DATA` and regenerates the available analyses into that data root.
It needs the existing composition/fourfold tables, structural contact map, cybrid scores and partial energy tables.
The run used a separate data root, selected by `MITOAGE_PIPELINES_DATA`, with linked read-only inputs and copied
output seeds; it did not overwrite the previous analyses. `rerun.log` contains the run output with machine-specific
paths replaced by `$MITOAGE_PIPELINES_DATA`. The original log remains in the remote review data root.
