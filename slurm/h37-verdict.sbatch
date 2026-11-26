#!/bin/bash
# H37 part A: the entry's bar applied to every composition test and rate, after h37-rates and h37-tests. CPU only.
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h37-verdict.sbatch
#SBATCH --job-name=mitoage-h37-verdict
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
export UV_NO_SYNC=1
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
uv run --frozen mitoage longevity rna-verdict
sha256sum data/output/longevity/rna/rna_verdict.tsv | sed 's/  /\t/'
