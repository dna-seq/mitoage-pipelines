# The manuscript's analyses on the build with the EURING birds

Every analysis whose numbers the MitoAge manuscript reports (H1 `aa-pgls` in four classes, H2's column screen
`pgls` in mammals and birds, H7 `rates`, H8 `motif-pgls`, H11 `generation-time` in both maturity modes, H12
`structure-pgls`), run twice on haumea with `slurm/manuscript-rerun.sbatch` from commit 2ee27a1 of this repository:

| Run | Slurm job | Build (`public/data/` of dna-seq/mitoage) | species.parquet | dated tree (`phylogeny.json` at the same commit) |
|---|---|---|---|---|
| before | 6489876 | 4b5e433, 2,335 species (equal to longevity-genie/mitoage c525adbf) | 794c421dfcdc | 4b5e433 |
| after | 6489877 | fb54cec, 2,344 species: 9 European birds from the EURING list (D15 of the site repository), all genes realigned | 4e4aab66ae53 | fb54cec |

The build with the EURING birds is not on Hugging Face yet, so both runs read `public/data/` of the site's
repository at the named commit, copied into the data root (`MITOAGE_PIPELINES_DATA`); the logs record the digests.

- `*.tsv`: the tables of the "after" run that the manuscript reads.
- `comparison.txt`: `slurm/compare-runs.py before after`, table by table and class by class.
- `slurm-6489876.log.txt`, `slurm-6489877.log.txt`: the two jobs' logs.

**The "before" run reproduces the published tables**: byte for byte `generation_time.{female,any}{,.classes}.tsv`
(`docs/analyses/h11/`), `motif_pgls.tsv` (`results/h12/`), `rates.tsv` and `rate_contrasts.tsv`; and all 136 rows
of `results/h12/structure_pgls.tsv` (the run adds the second round's `mass+own composition` rows). So every
difference below comes from the build.

**What the nine birds change.** Two of the nine enter the classes' tested sets (birds 379 to 381; maturity sets 281
to 282). The other effect is the realignment: every gene was aligned again with the birds in it, so alignment
columns move for every class.

- Unchanged, byte for byte: `aa_pgls` in mammals, fishes and reptiles; every mammal row of `generation_time`,
  `motif_pgls` and `structure_pgls`; the fish rows of `motif_pgls` and `structure_pgls`; the reptile rows of
  `motif_pgls`.
- H2, mammals: 5,262 tests instead of 5,259 (columns moved); the same two columns pass, ND5 200 (q = 0.0050) and
  COX1 419 (q = 0.049). Birds: 2,738 tests, none passes (smallest p 0.0007, was 0.001).
- H1, birds (381): nothing at q < 0.05 (smallest p 0.018).
- H7: the same 339 cherries (101 bird pairs), 655 tests, 24.2% saturated; four tests pass, the same four (all
  species ND6 -0.197, q = 0.0075; birds CYTB per million years -0.022, q = 0.0081; mammals ND5 and ND2 unchanged).
  Birds on the concatenated genes per million years: p = 0.071, q = 0.28 (95 pairs). In mammals one cherry (the two
  elephants, in ATP8 and the concatenation) loses a codon to the realignment, which moves two of 140 mammal tests in
  the fourth digit; fish and reptile contrasts change the same way.
- H8: 132 tests, none passes (smallest q 0.16, unchanged).
- H11, birds (female maturity): maturity now follows 7 features alone, 8 with mass and 8 with lifespan and mass (was 7, 6, 7; codon
  CUU and the GC skew of the whole record join); lifespan still follows only UAU, in the 115 birds with a birth weight.
- H12, birds (381): contact minus control -0.015 (95% -0.037 to +0.008, q = 0.20). **Reptiles (97), with no reptile
  added: contact minus control -0.111 (-0.211 to -0.010, p = 0.031, q = 0.046), and the contact share -0.114
  (q = 0.046), where the build before gave q = 0.052.** The species are the same; the contact columns' residues
  moved with the realignment. A result that crosses q = 0.05 because a bird joined the alignment is not one to rely on.
