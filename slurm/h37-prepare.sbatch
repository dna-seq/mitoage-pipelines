#!/bin/bash
# H37 part A: the three small inputs of H1 and H7 (composition.tsv from `longevity controls`, fourfold.tsv, the dated tree),
# made once per data root before anything reads them. CPU only. From the root of a fresh clone of main whose data/input
# holds the build, and with H37_ROOTS naming the data roots to prepare, space separated (the clone's own data/ is the
# default root and needs no entry; others are folders laid out like data/ with their input/ linked):
#   uv sync --frozen
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h37-prepare.sbatch
#SBATCH --job-name=mitoage-h37-prepare
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
export UV_NO_SYNC=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
for root in "" $H37_ROOTS; do
  if [ -n "$root" ]; then export MITOAGE_PIPELINES_DATA="$root"; else unset MITOAGE_PIPELINES_DATA; fi
  echo "== data root ${root:-data/ of the clone}"
  sha256sum "${root:-data}/input/mitoage/species.parquet" "${root:-data}/input/mitoage-site/phylogeny.json"
  uv run --frozen mitoage longevity dated-tree
  uv run --frozen mitoage longevity controls > /dev/null
  uv run --frozen mitoage longevity fourfold
  sha256sum "${root:-data}"/output/longevity/composition.tsv "${root:-data}"/output/longevity/fourfold.tsv "${root:-data}"/interim/tree/dated.nwk
done
