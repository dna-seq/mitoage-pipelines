#!/bin/bash
# H37 part A, step 1: the 22 tRNA and 2 rRNA genes of every live species' record, one alignment per tRNA gene (MAFFT, the
# build's rule), the census of what each record carries, the features. From a fresh clone of main whose data/input holds
# the build, the pinned RefSeq release 237 under build2/refseq/release237 and the CONTIG-only records under
# build2/genbank/contig. H37_TOOLS names the folder that holds MAFFT. CPU only.
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --export=ALL,H37_TOOLS=<folder> --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h37-genes.sbatch
#SBATCH --job-name=mitoage-h37-genes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${H37_TOOLS:?set H37_TOOLS to the folder that holds mafft}"
export PATH="$H37_TOOLS:$PATH" UV_NO_SYNC=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
mafft --version 2>&1 | head -1
uv run --frozen mitoage longevity rna-genes --workers "$SLURM_CPUS_PER_TASK"
O=data/output/longevity/rna
sha256sum "$O"/*.tsv "$O"/*.fasta | sed 's/  /\t/'
