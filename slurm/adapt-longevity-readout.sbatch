#!/bin/bash
# H26's readout, from a clone of main at the corrected pair code (4c3f9fc or later): the trait and the pairs again,
# checked against the committed corrected tables, then the finished RELAX and aBSREL runs of
# slurm/adapt-longevity-selection.sbatch taken from the clone they ran in, then `adapt longevity-readout`.
#   sbatch --account=... --export=ALL,SEL=<that clone>/data/output/adapt-longevity slurm/adapt-longevity-readout.sbatch
# R1 and R4 do not use the corrected distances: the job prepares their inputs again and refuses when any differs from
# what the selection runs read. C is not run: it waits for H23's set.
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
  for f in pairs.tsv pair_tests.tsv; do
    cmp $O/$c/$f results/review-20261118/adapt-longevity/$c/$f || { echo "$c/$f differs from the committed corrected table"; exit 1; }
  done
  $M adapt longevity-selection-prep --class $c --out $O
  S=$O/$c/selection; T="$SEL/$c/selection"
  # every input HyPhy read there is the same bytes as the one prepared here; only its outputs may be extra
  extra=$(diff -rq $S "$T" | grep -vE "^Only in $T[^ ]*: (relax\.json|absrel\.(json|nwk|tsv))$" || true)
  [ -z "$extra" ] || { echo "selection inputs of $c differ from the runs' clone:"; echo "$extra"; exit 1; }
  for d in "$T"/*/draw*/; do
    r=${d#"$T"/}
    for f in relax.json absrel.json absrel.nwk absrel.tsv; do [ -f "$d$f" ] && cp "$d$f" $S/$r$f; done
  done
done
echo "pairs: the committed corrected tables' bytes; selection inputs: the same bytes as the ones the runs read"
$M adapt longevity-readout --out $O
sha256sum $O/readout.tsv $O/*/pair_tests.tsv
