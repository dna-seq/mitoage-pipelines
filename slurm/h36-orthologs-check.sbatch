#!/bin/bash
# H36, donor orthologs from RefSeq: the tests of the swap, of its model and of the donor orthologs, then the breaks of each guard
# of the donor orthologs and of the swap's reading of a donor's protein (slurm/h36-orthologs-mutations.py) and of the swap
# (slurm/h36-swap-mutations.py). CPU only; no network. From the root of a clone of main that holds the inputs of the swap (the
# structures, the cybrid records, the build, the vertebrate alignments, the weights) and `mitoage download donor-orthologs`:
#   UV_PROJECT_ENVIRONMENT=.venv-stab uv sync --frozen --group stab --group dev
#   mkdir -p slurm/logs && sbatch --account=<account> --nice=100 --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h36-orthologs-check.sbatch
#SBATCH --job-name=mitoage-h36-orth-check
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
export UV_PROJECT_ENVIRONMENT=.venv-stab
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
R="uv run --frozen --no-sync"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
sha256sum data/input/donor-orthologs/*.faa src/mitoage_pipelines/mitoport/donor_orthologs.tsv
$R pytest -q -p no:cacheprovider -rs tests
set +e
$R python slurm/h36-orthologs-mutations.py; a=$?
$R python slurm/h36-swap-mutations.py; b=$?
echo "breaks: orthologs exit $a, swap exit $b"
[ $a -eq 0 ] && [ $b -eq 0 ]
