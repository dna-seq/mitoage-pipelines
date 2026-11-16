#!/bin/bash
# H18: the interface energy of each cybrid's chimeric complex (`mitoport cybrid-energy`), one array task per (host,
# donor, complex), on CPU nodes; no GPU is used. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main (network is needed there, not in the job):
#   UV_PROJECT_ENVIRONMENT=.venv-energy uv sync --frozen --group energy
#   UV_PROJECT_ENVIRONMENT=.venv-energy uv run --frozen --group energy mitoage download pdb 5XTD 5XTE 5Z62 8H9S 8OM1 7O3H 7O37
#   UV_PROJECT_ENVIRONMENT=.venv-energy uv run --frozen --group energy mitoage download cybrid
#   mkdir -p slurm/logs && sbatch slurm/h18-energy.sbatch      # tasks: `mitoage mitoport cybrid-energy --list`
#   sbatch --dependency=afterany:<array job> slurm/h18-energy-check.sbatch
#SBATCH --job-name=mitoage-h18-energy
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --array=0-43%12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-energy
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); task $SLURM_ARRAY_TASK_ID on $(hostname)"
sha256sum data/input/cybrid/records.gb data/input/pdb/{5XTD,5XTE,5Z62,8H9S,8OM1,7O3H,7O37}.cif
$UV run --frozen --group energy mitoage mitoport cybrid-energy --task "$SLURM_ARRAY_TASK_ID" --threads "$SLURM_CPUS_PER_TASK"
