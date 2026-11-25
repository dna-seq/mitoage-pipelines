#!/bin/bash
# H36, step 3, the probe: SaProt 650M loaded on one idle card of another group's node (the requeue partition; the
# A100 cards are never used for this work) and one real pair read (`mitoage mitoport saprot --probe`): which attention it
# runs with on a Turing card, the time and memory of a read, the PLL with the real and a shuffled partner, and whether a read repeated
# is the same. It writes nothing. From the root of a fresh clone of main that holds the inputs (the structures, the cybrid records, the
# weights of `mitoage download saprot` and the tokens of `mitoage mitoport saprot-tokens`):
#   UV_PROJECT_ENVIRONMENT=.venv-saprot uv sync --frozen --group saprot --group stab --group dev
#   mkdir -p slurm/logs && sbatch --comment 'seat=mitoage/lead; what=H36 step 3 probe; eta=?; cuda=?' slurm/h36-step3-probe.sbatch
#SBATCH --job-name=mitoage-h36-probe
#SBATCH --partition=GPUrequeue
#SBATCH --account=cluster
#SBATCH --constraint=gpu
#SBATCH --nodelist=gpu1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
source slurm/h36-step3-card.sh
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-saprot HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD) on $(hostname)"
nvidia-smi --query-gpu=index,name,memory.total,memory.used,utilization.gpu,compute_cap --format=csv
nvidia-smi --query-compute-apps=gpu_uuid,pid,used_memory --format=csv
take_card
nvidia-smi --query-gpu=index,name,memory.used --format=csv -i "$CUDA_VISIBLE_DEVICES"
$UV run --frozen --no-sync mitoage mitoport saprot --probe --device cuda
