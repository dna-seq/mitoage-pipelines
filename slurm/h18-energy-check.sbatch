#!/bin/bash
# H18, after the array of slurm/h18-energy.sbatch: collect the tasks' tables, H16's baseline and the readout, then the
# tests and the breaks of each guard, in the environment with OpenMM. CPU only.
#   sbatch --dependency=afterany:<array job> slurm/h18-energy-check.sbatch
#SBATCH --job-name=mitoage-h18-check
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-energy
R="$UV run --frozen --group energy --group dev"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
O=data/output/mitoport/cybrid
$R mitoage mitoport cybrid > /dev/null
$R mitoage mitoport cybrid-energy-collect
$R mitoage mitoport cybrid-energy-readout
sha256sum $O/energy_scores.tsv $O/h18_calibration.tsv $O/h18_summary.tsv
$R pytest -q -p no:cacheprovider tests/test_cybrid_energy.py tests/test_cybrid.py tests/test_cybrid_h17.py
$R python slurm/h18-energy-mutations.py
