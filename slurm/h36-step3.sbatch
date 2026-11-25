#!/bin/bash
# H36, step 3: SaProt 650M on the joined mitochondrial and nuclear chains (`mitoage mitoport saprot`), one array element per shard of the
# swap's tasks (every 16th task from the element's number), each task writing its own tables, so that a requeued element (the requeue
# partition's jobs are taken back without warning when the node's owner needs the card) starts again and skips the tasks that have their
# table. Idle cards of another group's node only, never the A100 cards; the card is taken under a node-local lock
# (slurm/h36-step3-card.sh). From the root of a fresh clone of main that holds the inputs (see slurm/h36-step3-probe.sbatch):
#   sbatch --comment 'seat=mitoage/lead; what=H36 step 3 SaProt human host; eta=?; cuda=?' --array=0-15%4 --export=ALL,H36_HOST=homo-sapiens,H36_OUT=<folder> slurm/h36-step3.sbatch
#   sbatch ... --export=ALL,H36_HOST=mus-musculus,H36_OUT=<folder> slurm/h36-step3.sbatch
#SBATCH --job-name=mitoage-h36-saprot
#SBATCH --partition=GPUrequeue
#SBATCH --account=cluster
#SBATCH --constraint=gpu
#SBATCH --nodelist=gpu1
#SBATCH --array=0-15%4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
source slurm/h36-step3-card.sh
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-saprot HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
HOST=${H36_HOST:-homo-sapiens}
OUT=${H36_OUT:-data/output/mitoport/cybrid-saprot-$HOST}
export H36_COMMENT_PREFIX="seat=mitoage/lead; what=H36 step 3 SaProt $HOST shard $SLURM_ARRAY_TASK_ID/16"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); shard $SLURM_ARRAY_TASK_ID of 16 on $(hostname); host $HOST; out $OUT"
sha256sum data/input/cybrid/records.gb data/input/pdb/{5XTD,5XTE,5Z62,8H9S,7O37,7O3H,8OM1}.cif data/input/mitoage/species.parquet data/interim/saprot-tokens/*.tsv 2>/dev/null || true
take_card
s=$(date +%s)
$UV run --frozen --no-sync mitoage mitoport saprot --host "$HOST" --shard "$SLURM_ARRAY_TASK_ID/16" --out "$OUT" --device cuda
echo "shard $SLURM_ARRAY_TASK_ID: $(( $(date +%s) - s )) s"
