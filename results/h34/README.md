# H34, part A: the predictors that already exist, on H29's variants

Run on the login node from a fresh clone at commit 832eca5 (a table join, seconds): `download mitimpact`
(MitImpact 4, archive and table pinned by SHA-256), then `mitoport variant-compare` on H29's corrected
`variant_scores.tsv` and H31's corrected `kron_variants.tsv` (`results/review-20261122/`). Tests 7 passed; breaks
10 of 10 caught. Reported, claimed as nothing: `docs/HYPOTHESES.md`, H34.

`h34_existing.tsv`: one row per score, with the variants it covers, its area under the ROC curve for harmful
against harmless, the area of the site's variability among the build's mammals on the same variants, the
difference with its 95% bootstrap interval, and what the score knew of labels like ours.

| score | harmful, harmless | area | baseline there | difference (interval) | labels |
| --- | --- | --- | --- | --- | --- |
| APOGEE 2 | 39, 328 | 1.000 | 0.919 | 0.081 (0.030 to 0.146) | trained on such variants of these genes: not a test |
| AlphaMissense | 39, 325 | 0.979 | 0.918 | 0.061 (0.015 to 0.117) | calibrated with population frequency |
| ESM1b | 39, 313 | 0.937 | 0.919 | 0.018 (-0.057 to 0.095) | none |
| KronGenome (H31) | 39, 328 | 0.923 | 0.919 | 0.004 (-0.058 to 0.072) | none |
| PhyloP, 100 vertebrates | 39, 328 | 0.838 | 0.919 | -0.081 (-0.155 to -0.010) | none |
| gnomAD local constraint | 39, 328 | 0.735 | 0.919 | -0.184 (-0.306 to -0.071) | built from population variation |
| StaB-ddG (H29) | 38, 314 | 0.629 | 0.915 | -0.286 (-0.422 to -0.150) | none |

Reading: the two scores that saw no labels and read sequence alone (ESM1b, KronGenome) are level with counting
the amino acids a site holds among mammals. Only the two that were fitted to labels like ours are above it, and
that is what such fitting gives. The bar a new model has to clear is therefore the baseline's 0.92, as H34 says.
