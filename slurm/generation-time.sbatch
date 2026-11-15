#!/bin/bash
# H11 (`mitoage longevity generation-time`) on the cluster host haumea, with its tests and a rerun check.
# The cluster runs only code that is on main. On the login node, from the root of a checkout pulled from main:
#   uv sync --frozen
#   uv run mitoage download build --revision c525adbff8c12b19981fdada7f1615c2559b7140
#   uv run mitoage download tree && uv run mitoage longevity dated-tree && uv run mitoage longevity fourfold
#   mkdir -p slurm/logs && sbatch slurm/generation-time.sbatch
# Writes data/output/longevity/generation_time.{female,any}{,.classes}.tsv and slurm/logs/<job>.log, which
# records the commit, the tests, the two runs and whether a second run reproduces the first byte for byte.
#SBATCH --job-name=mitoage-h11
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
# git is not installed on the compute nodes: read the commit from .git directly.
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet"
$UV run --frozen pytest -q tests/test_generation_time.py
out=data/output/longevity
for m in female any; do
  $UV run --frozen mitoage longevity generation-time --maturity "$m"
  cp "$out/generation_time.$m.tsv" "$out/generation_time.$m.first.tsv"
  $UV run --frozen mitoage longevity generation-time --maturity "$m" > /dev/null
  if cmp -s "$out/generation_time.$m.tsv" "$out/generation_time.$m.first.tsv"; then
    echo "rerun $m: identical bytes"
  else
    echo "rerun $m: DIFFERS"
  fi
  rm "$out/generation_time.$m.first.tsv"
done
sha256sum "$out"/generation_time.*.tsv
