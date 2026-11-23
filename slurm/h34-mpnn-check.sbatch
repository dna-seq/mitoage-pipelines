#!/bin/bash
# H34, part B, after the array of slurm/h34-mpnn.sbatch: collect the complexes' tables with the baselines, the
# readout, the byte comparison with a second run when H34_SECOND names its folder, the tests, the breaks, and the
# check of our reading against the LigandMPNN repository's own score.py (which needs ProDy, added for that step only).
# CPU only.
#   sbatch --account=<account> --nice=100 --dependency=afterany:<array job> [--export=ALL,H34_SECOND=<folder>] slurm/h34-mpnn-check.sbatch
#SBATCH --job-name=mitoage-h34-check
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
O=data/output/mitoport/variants-mpnn
$R mitoage mitoport variant-mpnn-collect
$R mitoage mitoport variant-mpnn-readout
sha256sum $O/mpnn_scores.tsv $O/h34_mpnn_summary.tsv
if [ -n "$H34_SECOND" ]; then
  same=0; differ=0; absent=0
  for f in $O/mpnn/*.tsv; do
    g="$H34_SECOND/mpnn/$(basename $f)"
    if [ ! -f "$g" ]; then absent=$((absent+1)); elif cmp -s "$f" "$g"; then same=$((same+1)); else differ=$((differ+1)); echo "differs: $(basename $f)"; fi
  done
  echo "second run: $same tables identical byte for byte, $differ differ, $absent not run a second time"
else
  echo "second run: not compared (no folder given)"
fi
$R pytest -q -p no:cacheprovider tests/test_variant_mpnn.py tests/test_variant_stab.py tests/test_cybrid_stab.py
$R python slurm/h34-mpnn-mutations.py
mkdir -p data/interim/h34-tool-check
$R --with prody python slurm/h34-mpnn-tool-check.py data/interim/h34-tool-check
