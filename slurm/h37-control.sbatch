#!/bin/bash
# H37 part A, the positive control: the rate pipeline of the RNA genes on the 13 protein genes against H7's table, row by
# row. Two tasks. Task 0: the build of H7's published numbers (0f4e970), in the data root H37_PUBLISHED_ROOT, against the
# repository's results/review-20261118/longevity/rates.tsv. Task 1: the build H37 runs on, in the clone's own data/, against
# H37_REFERENCE, H7's rates.tsv from the run of H27 on that build. Both after slurm/h37-prepare.sbatch. CPU only.
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --array=0-1 --export=ALL,H37_PUBLISHED_ROOT=<folder>,H37_REFERENCE=<rates.tsv> --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h37-control.sbatch
#SBATCH --job-name=mitoage-h37-control
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${H37_PUBLISHED_ROOT:?}" "${H37_REFERENCE:?}"
export UV_NO_SYNC=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
case "$SLURM_ARRAY_TASK_ID" in
  0) export MITOAGE_PIPELINES_DATA="$H37_PUBLISHED_ROOT"; REF=results/review-20261118/longevity/rates.tsv; D="$H37_PUBLISHED_ROOT" ;;
  1) unset MITOAGE_PIPELINES_DATA; REF="$H37_REFERENCE"; D=data ;;
  *) echo "no task $SLURM_ARRAY_TASK_ID"; exit 2 ;;
esac
sha256sum "$D/input/mitoage/species.parquet" "$REF"
uv run --frozen mitoage longevity rna-control --reference "$REF"
sha256sum "$D/output/longevity/rna/h7_control_rates.tsv" "$D/output/longevity/rna/h7_control.json"
