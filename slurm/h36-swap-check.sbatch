#!/bin/bash
# H36, step 2, after the array of slurm/h36-swap.sbatch: collect the tasks' tables into the swap tables and apply the
# candidate rule, the byte comparison with a second run when H36_SECOND names its folder, the tests and the breaks of
# each guard. CPU only. The environment is not synced here (--no-sync): other tasks may still be using it.
#   sbatch --account=<account> --nice=100 --comment '...' --dependency=afterany:<array job> [--export=ALL,H36_SECOND=<folder>] slurm/h36-swap-check.sbatch
#SBATCH --job-name=mitoage-h36-check
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
R="$UV run --frozen --no-sync"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
O=data/output/mitoport/cybrid
$R mitoage mitoport swap-collect
sha256sum $O/swap_*.tsv
if [ -n "$H36_SECOND" ]; then
  $R mitoage mitoport swap-collect --out "$H36_SECOND" > /dev/null
  same=0; differ=0; absent=0
  for f in $O/swap/*.tsv; do
    g="$H36_SECOND/swap/$(basename $f)"
    if [ ! -f "$g" ]; then absent=$((absent+1)); elif cmp -s "$f" "$g"; then same=$((same+1)); else differ=$((differ+1)); echo "differs: $(basename $f)"; fi
  done
  echo "second run: $same task tables identical byte for byte, $differ differ, $absent not run a second time"
  for f in $O/swap_*.tsv; do
    g="$H36_SECOND/$(basename $f)"
    if cmp -s "$f" "$g"; then echo "identical: $(basename $f)"; else echo "differs: $(basename $f)"; fi
  done
else
  echo "second run: not compared (no folder given)"
fi
$R pytest -q -p no:cacheprovider tests/test_cybrid_swap.py tests/test_cybrid_stab.py
$R python slurm/h36-swap-mutations.py
