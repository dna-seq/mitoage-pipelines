#!/bin/bash
# H28 and H29 on one GPU card, when one is free: every task of `mitoport cybrid-stab` and, with H29=1, of `mitoport
# variant-stab`, one after another, each printing its time and its peak GPU memory. A task whose table is there is
# not run again. The card is named by the submitter, who has looked that it is free (no job is put on a card in use):
#   sbatch --account=<account> --partition=<GPU partition> --nodelist=<node> --export=ALL,CARD=<card>[,H29=1][,H28_OUT=<folder>][,H29_OUT=<folder>] slurm/h28-stab-gpu.sbatch
# Several cards: one job per card, each with SHARD=<k>,SHARDS=<n> (k from 0 to n-1); a job runs the tasks whose
# number leaves k when divided by n. H28=0 leaves H28's tasks out.
#SBATCH --job-name=mitoage-h28-stab-gpu
#SBATCH --constraint=A100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
if [ -z "$CARD" ]; then echo "CARD is not set: name the free card"; exit 2; fi
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-stab CUDA_VISIBLE_DEVICES=$CARD
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
R="$UV run --frozen --group stab"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD) on $(hostname), card $CARD"
nvidia-smi --query-gpu=index,name,memory.used,memory.free --format=csv
sha256sum data/input/cybrid/records.gb data/input/pdb/*.cif data/input/stabddg/*.pt
K=${SHARD:-0}; N=${SHARDS:-1}
echo "shard $K of $N"
O=${H28_OUT:-data/output/mitoport/cybrid}
[ "${H28:-1}" = 0 ] || $R mitoage mitoport cybrid-stab --list | while read t host donor cx; do
  [ $((t % N)) -eq "$K" ] || continue
  if [ -f "$O/stab/$host.$donor.$cx.tsv" ]; then echo "task $t: table there, not run again"; continue; fi
  s=$(date +%s)
  $R mitoage mitoport cybrid-stab --task "$t" --device cuda --out "$O" < /dev/null
  echo "task $t ($host $donor $cx): $(( $(date +%s) - s )) s"
done
if [ -n "$H29" ]; then
  sha256sum data/input/variants/*
  V=${H29_OUT:-data/output/mitoport/variants}
  for t in $(seq 0 18); do
    [ $((t % N)) -eq "$K" ] || continue
    s=$(date +%s)
    $R mitoage mitoport variant-stab --task "$t" --device cuda --out "$V" < /dev/null
    echo "variants task $t: $(( $(date +%s) - s )) s"
  done
fi
