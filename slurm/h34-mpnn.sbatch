#!/bin/bash
# H34, part B: ProteinMPNN's membrane checkpoint and its three controls on H29's human missense variants
# (`mitoport variant-mpnn`), one array task per complex, on CPU nodes (the reading is light; a GPU is not needed).
# The cluster runs only code that is on main. On the login node, from the root of a fresh clone of main that holds
# H29's inputs (data/input/cybrid/records.gb, data/input/pdb, data/input/variants as fetched for H29, the build in
# data/input/mitoage, the sites of `mitoport burial` in data/interim/mitoport) and, in data/input/ligandmpnn, the
# LigandMPNN repository at the pinned commit with its weights (`get_model_params.sh`):
#   UV_PROJECT_ENVIRONMENT=.venv-stab uv sync --frozen --group stab --group dev
#   mkdir -p slurm/logs && sbatch --account=<account> --nice=100 slurm/h34-mpnn.sbatch
#   sbatch --account=<account> --nice=100 --dependency=afterany:<array job> [--export=ALL,H34_SECOND=<folder>] slurm/h34-mpnn-check.sbatch
# Every task a second time into another folder, for the byte comparison of the check job:
#   sbatch --account=<account> --nice=100 --export=ALL,H34_OUT=<folder> slurm/h34-mpnn.sbatch
#SBATCH --job-name=mitoage-h34-mpnn
#SBATCH --partition=compute
#SBATCH --array=0-3
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
sha256sum data/input/cybrid/records.gb data/input/pdb/{5XTD,5XTE,5Z62,8H9S}.cif data/input/variants/* data/input/ligandmpnn/model_utils.py data/input/ligandmpnn/model_params/*.pt
s=$(date +%s)
# --no-sync: the tasks start together on one environment, and a `uv run` that syncs it (to other groups than the
# ones installed) removes and rewrites torch's files under the others; four of eight tasks died of it on 2026-11-23
$UV run --frozen --no-sync mitoage mitoport variant-mpnn --task "$SLURM_ARRAY_TASK_ID" ${H34_OUT:+--out "$H34_OUT"}
echo "task $SLURM_ARRAY_TASK_ID: $(( $(date +%s) - s )) s"
