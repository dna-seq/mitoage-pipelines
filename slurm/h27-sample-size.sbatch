#!/bin/bash
# H27, sample size as a covariate (docs/HYPOTHESES.md): the lifespan slopes of H1 and H11 without and with AnAge's
# sample size in the model, on the set not_unreliable. It reads the data root of that set as slurm/h27-sets.sbatch
# left it (the dated tree, the composition and fourfold tables) and writes into <base>/sample_size.
#   sbatch --account=<account> --export=ALL,H27_BASE=<base> slurm/h27-sample-size.sbatch
#SBATCH --job-name=mitoage-h27-sample-size
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${H27_BASE:?set H27_BASE to the folder that holds the data roots of the H27 run}"
export MITOAGE_PIPELINES_DATA="$H27_BASE/not_unreliable"
export MITOAGE_SPECIES_SET=not_unreliable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
UV=${UV:-uv}
# git is not installed on the compute nodes: read the commit from .git directly.
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
sha256sum "$MITOAGE_PIPELINES_DATA/input/mitoage/species.parquet" "$MITOAGE_PIPELINES_DATA/output/longevity/composition.tsv"
$UV run --frozen mitoage longevity sample-size-covariate --out "$H27_BASE/sample_size"
