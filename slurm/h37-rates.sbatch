#!/bin/bash
# H37 part A: the RNA genes' rates against lifespan (H7's relative-rate test with the distance of a stretch with no codons),
# and the 13 protein genes' rates for the bar's comparison, on H7's species and on H27's high set. Refused by the code
# unless the control of this build reproduced H7. After h37-genes and h37-control. CPU only.
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h37-rates.sbatch
#SBATCH --job-name=mitoage-h37-rates
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
export UV_NO_SYNC=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
uv run --frozen mitoage longevity rna-rates
O=data/output/longevity/rna
sha256sum "$O"/rna_rates.tsv "$O"/protein_rates.tsv "$O"/rna_rate_contrasts.*.tsv | sed 's/  /\t/'
