#!/bin/bash
# H17, GPU part: ESM-2 650M masked marginals of the donors' residues in the hosts' proteins, alone and beside the
# host's nuclear partner (`mitoport cybrid-esm`), then the readout of every H17 score. Submit after the CPU part:
#   sbatch --dependency=afterok:<h17-cybrid job> slurm/h17-cybrid-esm.sbatch
# On the login node first (network is needed there, not in the job), in its own environment as for H12:
#   UV_PROJECT_ENVIRONMENT=.venv-fold uv sync --frozen --group fold
#   UV_PROJECT_ENVIRONMENT=.venv-fold HF_HOME=/data/uq2347/.cache/huggingface uv run --frozen --group fold python -c \
#     "from transformers import AutoModelForMaskedLM as A, AutoTokenizer as T; m='facebook/esm2_t33_650M_UR50D'; T.from_pretrained(m); A.from_pretrained(m)"
#SBATCH --job-name=mitoage-h17-esm
#SBATCH --partition=gpu1024
#SBATCH --account=itmz
#SBATCH --constraint=A100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
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
O=data/output/mitoport/cybrid
$R mitoage mitoport cybrid-esm
cat $O/esm_checkpoint.txt
$R mitoage mitoport cybrid-readout
sha256sum $O/*.tsv
