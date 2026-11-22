#!/bin/bash
# H29: StaB-ddG's predicted change of stability for human missense variants of the mtDNA-encoded proteins, harmful
# against harmless (`mitoport variant-stab`), one array task per complex and share of its variants, on CPU nodes.
# The cluster runs only code that is on main. On the login node, from the root of a fresh clone of main that holds
# H28's inputs (the structures, the records, the sites, the weights) and the build (network is needed there only):
#   UV_PROJECT_ENVIRONMENT=.venv-stab uv sync --frozen --group stab --group dev
#   UV_PROJECT_ENVIRONMENT=.venv-stab uv run --frozen --group stab mitoage download variants
#   UV_PROJECT_ENVIRONMENT=.venv-stab uv run --frozen --group stab mitoage mitoport variant-stab --list
#   mkdir -p slurm/logs && sbatch --account=<account> slurm/h29-variants.sbatch
#   sbatch --account=<account> --dependency=afterany:<array job> slurm/h29-variants-check.sbatch
# Every task a second time into another folder, for the byte comparison of the check job:
#   sbatch --account=<account> --export=ALL,H29_OUT=<folder> slurm/h29-variants.sbatch
#SBATCH --job-name=mitoage-h29-variants
#SBATCH --partition=compute
#SBATCH --array=0-18%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-stab
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); task $SLURM_ARRAY_TASK_ID on $(hostname)"
sha256sum data/input/cybrid/records.gb data/input/pdb/{5XTD,5XTE,5Z62,8H9S}.cif data/input/stabddg/*.pt data/input/variants/*
$UV run --frozen --group stab mitoage mitoport variant-stab --task "$SLURM_ARRAY_TASK_ID" ${H29_OUT:+--out "$H29_OUT"}
