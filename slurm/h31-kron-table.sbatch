#!/bin/bash
# H31: KronGenome's likelihood tables of the two hosts' mitochondrial genomes (three arms each), twice, the two sets
# compared byte for byte. On a GPU card name it with CARD; without CARD it runs on the CPU. WEIGHTS is the folder with
# config.json and model.safetensors of the pinned revision.
#   sbatch --account=<account> --partition=<partition> --export=ALL,WEIGHTS=<folder>[,CARD=<card>] slurm/h31-kron-table.sbatch
#SBATCH --job-name=mitoage-h31-kron
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
if [ -z "$WEIGHTS" ]; then echo "WEIGHTS is not set"; exit 2; fi
UV=${UV:-uv}
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
D=cpu; if [ -n "$CARD" ]; then export CUDA_VISIBLE_DEVICES=$CARD; D=cuda; fi
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD) on $(hostname), device $D"
O=${H31_OUT:-data/output/mitoport/kron}
for run in tables tables-second; do
  $UV run --script slurm/h31-kron-table.py --weights "$WEIGHTS" --gb data/input/cybrid/records.gb --out "$O/$run" --device $D < /dev/null
done
for f in "$O"/tables/*.tsv; do cmp "$f" "$O/tables-second/$(basename "$f")" && echo "same bytes: $(basename "$f")"; done
sha256sum "$O"/tables/*.tsv
