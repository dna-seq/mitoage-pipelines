#!/bin/bash
# H28, after the array of slurm/h28-stab.sbatch: collect the tasks' tables, H16's baseline and the readout, compare
# the tables with those of a second run when H28_SECOND names its folder, then the tests and the breaks of each guard,
# in the environment with the model. CPU only.
#   sbatch --dependency=afterany:<array job>[:<second array job>] [--export=ALL,H28_SECOND=<folder>] slurm/h28-stab-check.sbatch
#SBATCH --job-name=mitoage-h28-check
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
O=data/output/mitoport/cybrid
$R mitoage mitoport cybrid > /dev/null
$R mitoage mitoport cybrid-stab-collect
$R mitoage mitoport cybrid-stab-readout
sha256sum $O/stab_scores.tsv $O/stab_by_seed.tsv $O/h28_calibration.tsv $O/h28_summary.tsv $O/h28_verdict.tsv
if [ -n "$H28_SECOND" ]; then
  same=0; differ=0; absent=0
  for f in $O/stab/*.tsv; do
    g="$H28_SECOND/stab/$(basename $f)"
    if [ ! -f "$g" ]; then absent=$((absent+1)); elif cmp -s "$f" "$g"; then same=$((same+1)); else differ=$((differ+1)); echo "differs: $(basename $f)"; fi
  done
  echo "second run: $same tables identical byte for byte, $differ differ, $absent not run a second time"
else
  echo "second run: not compared (no folder given)"
fi
$R pytest -q -p no:cacheprovider tests/test_cybrid_stab.py tests/test_cybrid.py tests/test_cybrid_h17.py tests/test_cybrid_energy.py
$R python slurm/h28-stab-mutations.py
