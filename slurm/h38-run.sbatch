#!/bin/bash
# H38: the architecture of the mitochondrial genome against lifespan, from a fresh clone of main whose data/input holds the build (with
# data_confidence), the dated species tree's source, RefSeq release 237 under build2/refseq/release237, the CONTIG-only records under
# build2/genbank/contig and NCBI Taxonomy under build2/taxdump/current-2026-11-13. H38_TOOLS names the folder that holds Tandem Repeats
# Finder 4.09.1 (bioconda build h031d066_4; any other binary is refused). In order: the audit's records, strand, order and code (H13), the
# small inputs of H1 (composition.tsv and the dated tree), the measures, the positive control (the job stops there, with a non-zero exit,
# if the rearrangement measure does not reproduce the audit's flags), the tests. CPU only.
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --export=ALL,H38_TOOLS=<folder> --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h38-run.sbatch
#SBATCH --job-name=mitoage-h38-run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${H38_TOOLS:?set H38_TOOLS to the folder that holds trf}"
export PATH="$H38_TOOLS:$PATH" UV_NO_SYNC=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
sha256sum data/input/mitoage/species.parquet data/input/mitoage-site/phylogeny.json "$H38_TOOLS/trf"
M="uv run --frozen mitoage"
$M audit records
$M audit strand
$M audit order
$M audit code
$M longevity dated-tree
$M longevity controls > /dev/null
$M longevity arch-measures --workers "$SLURM_CPUS_PER_TASK"
$M longevity arch-control --published results/job3-audit/0f4e970/order_flags.tsv
$M longevity arch-tests
O=data/output/longevity/arch
A=data/output/audit
sha256sum "$O"/*.tsv "$O"/*.json "$A"/order.tsv "$A"/order_flags.tsv "$A"/strand.tsv "$A"/code.tsv "$A"/annotation_vs_build.tsv data/output/longevity/composition.tsv data/interim/tree/dated.nwk | sed 's/  /\t/'
