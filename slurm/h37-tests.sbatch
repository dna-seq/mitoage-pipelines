#!/bin/bash
# H37 part A: the RNA genes' composition against lifespan in one taxon group per array task (0 Mammalia, 1 Aves, 2
# Actinopterygii, 3 Reptilia), with the 13 protein genes as the comparison. Refused by the code unless the control of this
# build reproduced H7. After h37-genes and h37-control. CPU only.
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --array=0-3 --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h37-tests.sbatch
#SBATCH --job-name=mitoage-h37-tests
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
export UV_NO_SYNC=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
CLASSES=(Mammalia Aves Actinopterygii Reptilia)
C=${CLASSES[$SLURM_ARRAY_TASK_ID]:?no class for array task $SLURM_ARRAY_TASK_ID}
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); class $C"
uv run --frozen mitoage longevity rna-tests --class "$C"
sha256sum data/output/longevity/rna/rna_tests."$C".tsv | sed 's/  /\t/'
