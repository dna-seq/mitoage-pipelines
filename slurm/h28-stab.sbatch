#!/bin/bash
# H28: StaB-ddG's predicted change of binding energy for each cybrid's complex (`mitoport cybrid-stab`), one array task
# per (host, donor, complex), on CPU nodes; no GPU is used. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main that holds H17's inputs (the structures, the records and
# the sites of `mitoport burial`; network is needed there, not in the job):
#   UV_PROJECT_ENVIRONMENT=.venv-stab uv sync --frozen --group stab --group dev
#   UV_PROJECT_ENVIRONMENT=.venv-stab uv run --frozen --group stab mitoage download stabddg
#   mkdir -p slurm/logs && sbatch --account=<account> slurm/h28-stab.sbatch       # tasks: `mitoage mitoport cybrid-stab --list`
#   sbatch --dependency=afterany:<array job> slurm/h28-stab-check.sbatch
# Every task a second time into another folder, for the byte comparison of the check job:
#   sbatch --export=ALL,H28_OUT=<folder> slurm/h28-stab.sbatch
#SBATCH --job-name=mitoage-h28-stab
#SBATCH --partition=compute
#SBATCH --array=0-43%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-stab
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); task $SLURM_ARRAY_TASK_ID on $(hostname)"
sha256sum data/input/cybrid/records.gb data/input/pdb/{5XTD,5XTE,5Z62,8H9S,8OM1,7O3H,7O37}.cif data/input/stabddg/*.pt
$UV run --frozen --group stab mitoage mitoport cybrid-stab --task "$SLURM_ARRAY_TASK_ID" ${H28_OUT:+--out "$H28_OUT"}
