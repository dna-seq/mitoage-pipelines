#!/bin/bash
# H26's readout, from a clone of main at the corrected pair code (4c3f9fc or later): the trait and the pairs again,
# checked against the committed corrected tables, then the finished RELAX and aBSREL runs of
# slurm/adapt-longevity-selection.sbatch taken from the clone they ran in, then `adapt longevity-readout`.
#   sbatch --account=... --export=ALL,SEL=<that clone>/data/output/adapt-longevity slurm/adapt-longevity-readout.sbatch
# R1 and R4 take their species from the pair table and not from the corrected distances; the job refuses when the
# pair table here differs from the one the selection runs were prepared from. C is not run: it waits for H23's set.
#SBATCH --job-name=mitoage-adapt-longevity-readout
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 POLARS_MAX_THREADS=4
[ -d "$SEL" ] || { echo "SEL is not a folder: set it to the adapt-longevity folder the selection runs wrote"; exit 1; }
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet"
$M longevity dated-tree
O=data/output/adapt-longevity
$M adapt longevity-traits --out $O
$M adapt longevity-pairs --out $O
for c in Mammalia Aves; do
  cmp $O/$c/pairs.tsv "$SEL/$c/pairs.tsv" || { echo "pairs of $c differ from the selection runs' clone"; exit 1; }
  for f in pairs.tsv pair_tests.tsv; do
    cmp $O/$c/$f results/review-20261118/adapt-longevity/$c/$f || { echo "$c/$f differs from the committed corrected table"; exit 1; }
  done
  rm -rf $O/$c/selection
  cp -r "$SEL/$c/selection" $O/$c/selection
done
echo "pairs: the same bytes as the selection runs' and as the committed corrected tables"
$M adapt longevity-readout --out $O
sha256sum $O/readout.tsv $O/*/pair_tests.tsv
