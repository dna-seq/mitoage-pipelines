#!/bin/bash
# H35, part 1: Evo 2's likelihood tables of the two hosts' mitochondrial genomes, twice, the two sets compared.
# Name the card with CARD and the checkpoint file with WEIGHTS:
#   sbatch --account=<account> --partition=<partition> --nodelist=<node> --export=ALL,WEIGHTS=<evo2_7b.pt>,CARD=<card> slurm/h35-evo2-table.sbatch
#SBATCH --job-name=mitoage-h35-evo2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
if [ -z "$WEIGHTS" ] || [ -z "$CARD" ]; then echo "WEIGHTS and CARD must be set"; exit 2; fi
UV=${UV:-uv}
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CUDA_VISIBLE_DEVICES=$CARD
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD) on $(hostname), card $CARD"
nvidia-smi --query-gpu=index,name,memory.used,memory.free --format=csv
O=${H35_OUT:-data/output/mitoport/evo2}
for run in tables tables-second; do
  $UV run --script slurm/h35-evo2-table.py --weights "$WEIGHTS" --gb data/input/cybrid/records.gb --out "$O/$run" < /dev/null
done
for f in "$O"/tables/*.tsv; do
  if cmp -s "$f" "$O/tables-second/$(basename "$f")"; then echo "same bytes: $(basename "$f")"
  else echo "differs: $(basename "$f"), largest difference $(paste "$f" "$O/tables-second/$(basename "$f")" | awk -F'\t' 'NR>1{for(i=3;i<=6;i++){d=$i-$(i+6); if(d<0)d=-d; if(d>m)m=d}} END{print m}')"; fi
done
sha256sum "$O"/tables/*.tsv
