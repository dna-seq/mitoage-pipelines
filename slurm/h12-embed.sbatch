#!/bin/bash
# H12, second round, GPU part: ESM-2 650M means of COX1-3 over complex IV's contact and control residues
# (`structure embed`) and their divergence against lifespan (`longevity structure-pgls --family embedding`), with the
# tests and the mutation check. Needs the inputs h12-partners.sbatch makes, so submit it after that job:
#   sbatch --dependency=afterok:<partners job> slurm/h12-embed.sbatch
# On the login node first (network is needed there, not in the job), in its own environment so that the CPU job's
# `uv run` does not remove torch from under it:
#   UV_PROJECT_ENVIRONMENT=.venv-fold uv sync --frozen --group fold
#   UV_PROJECT_ENVIRONMENT=.venv-fold HF_HOME=/data/uq2347/.cache/huggingface uv run --frozen --group fold python -c \
#     "from transformers import AutoModel, AutoTokenizer as T; m='facebook/esm2_t33_650M_UR50D'; T.from_pretrained(m); AutoModel.from_pretrained(m)"
#SBATCH --job-name=mitoage-h12-embed
#SBATCH --partition=gpu1024
#SBATCH --account=itmz
#SBATCH --constraint=A100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-fold HF_HOME=${HF_HOME:-/data/uq2347/.cache/huggingface} HF_HUB_OFFLINE=1
R="$UV run --frozen --group fold"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
nvidia-smi --query-gpu=index,name,memory.used,memory.free --format=csv
# never share a card that is busy: take one with at least 20 GiB free, else wait up to 2 h, then give up
for i in $(seq 1 24); do
  CARD=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits | sort -t, -k2 -rn | head -1)
  IDX=${CARD%%,*}; FREE=$(echo $CARD | cut -d, -f2 | tr -d ' ')
  if [ "$FREE" -ge 20000 ]; then echo "card $IDX has $FREE MiB free, using it"; export CUDA_VISIBLE_DEVICES=$IDX; break; fi
  echo "waiting: the freest card ($IDX) has only $FREE MiB free"; sleep 300
done
if [ -z "$CUDA_VISIBLE_DEVICES" ]; then echo "no free card on $(hostname) after 2 h, stopping"; exit 1; fi
$R mitoage structure embed
sha256sum data/interim/structure/embed_cox.parquet
$R pytest -q -p no:cacheprovider tests/test_structure_pgls.py
$R python slurm/structure-pgls-mutations.py embed
out=data/output/longevity/embedding
$R mitoage longevity structure-pgls --family embedding --out $out
cp $out/structure_pgls.tsv $out/structure_pgls.first.tsv
$R mitoage longevity structure-pgls --family embedding --out $out > /dev/null
if cmp -s $out/structure_pgls.tsv $out/structure_pgls.first.tsv; then echo "rerun: identical bytes"; else echo "rerun: DIFFERS"; fi
rm $out/structure_pgls.first.tsv
sha256sum $out/structure_pgls.tsv
