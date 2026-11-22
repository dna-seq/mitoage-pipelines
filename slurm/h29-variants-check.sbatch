#!/bin/bash
# H29, after the array of slurm/h29-variants.sbatch: collect the tasks' tables with the baselines, the readout, the
# byte comparison with a second run when H29_SECOND names its folder, then the tests and the breaks. CPU only.
#   sbatch --account=<account> --dependency=afterany:<array job> [--export=ALL,H29_SECOND=<folder>] slurm/h29-variants-check.sbatch
#SBATCH --job-name=mitoage-h29-check
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-stab
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
R="$UV run --frozen --group stab --group dev"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
O=data/output/mitoport/variants
$R mitoage mitoport variant-stab-collect
$R mitoage mitoport variant-stab-readout
sha256sum $O/variant_scores.tsv $O/h29_summary.tsv
if [ -n "$H29_SECOND" ]; then
  same=0; differ=0; absent=0
  for f in $O/variants/*.tsv; do
    g="$H29_SECOND/variants/$(basename $f)"
    if [ ! -f "$g" ]; then absent=$((absent+1)); elif cmp -s "$f" "$g"; then same=$((same+1)); else differ=$((differ+1)); echo "differs: $(basename $f)"; fi
  done
  echo "second run: $same tables identical byte for byte, $differ differ, $absent not run a second time"
else
  echo "second run: not compared (no folder given)"
fi
$R pytest -q -p no:cacheprovider tests/test_variant_stab.py tests/test_cybrid_stab.py
$R python slurm/h29-variants-mutations.py
