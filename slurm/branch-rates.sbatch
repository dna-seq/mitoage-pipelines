#!/bin/bash
# H3's two-step check of Coevol (`longevity branch-rates`, `branch-pgls`) on the cluster host haumea, run twice
# into two folders to check that a rerun gives the same tables. Needs `coevol-prep`'s inputs and the dated tree
# under $MITOAGE_PIPELINES_DATA, and hyphy on PATH (HyPhy 2.5.104).
#   mkdir -p slurm/logs && sbatch --export=ALL,MITOAGE_PIPELINES_DATA=<root> slurm/branch-rates.sbatch
#SBATCH --job-name=mitoage-branch-rates
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/branch-rates-%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${MITOAGE_PIPELINES_DATA:?}"
R="$MITOAGE_PIPELINES_DATA"; M="${UV:-uv} run --frozen mitoage"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); job $SLURM_JOB_ID on $(hostname)"
hyphy --version
sha256sum "$R"/interim/coevol/mammals.* "$R"/interim/tree/dated.nwk
for run in branch_rates branch_rates_rerun; do
  echo "== $run $(date -u +%FT%TZ)"
  /usr/bin/time -f "branch-rates %e s" $M longevity branch-rates --cpus "$SLURM_CPUS_PER_TASK" --out "$R/output/longevity/$run"
  $M longevity branch-pgls --out "$R/output/longevity/$run"
done
echo "== rerun identical?"
for f in branches.tsv model.json pgls.tsv; do
  cmp -s "$R/output/longevity/branch_rates/$f" "$R/output/longevity/branch_rates_rerun/$f" && echo "$f identical" || echo "$f DIFFERS"
done
sha256sum "$R"/output/longevity/branch_rates/*
