#!/bin/bash
# H37 part A: the tests of the RNA genes' code and the breaks of its guards (slurm/h37-mutations.py), then the rest of the
# suite, from a fresh clone of main after slurm/h37-prepare.sbatch. H37_TOOLS names the folder that holds MAFFT (the
# build's aligner, v7.526). CPU only.
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --export=ALL,H37_TOOLS=<folder> --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h37-check.sbatch
#SBATCH --job-name=mitoage-h37-check
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${H37_TOOLS:?set H37_TOOLS to the folder that holds mafft}"
export PATH="$H37_TOOLS:$PATH" UV_NO_SYNC=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
mafft --version 2>&1 | head -1
R="uv run --frozen"
set +e
$R pytest -q -p no:cacheprovider -rs tests/test_rna_genes.py tests/test_rna_rates.py tests/test_rna_tests.py; t=$?
echo "the tests of H37 part A: exit $t (a skipped test is listed above and counts as a failure of this check)"
$R pytest -q -p no:cacheprovider -rs tests/test_rna_genes.py tests/test_rna_rates.py tests/test_rna_tests.py 2>&1 | grep -c "^SKIPPED" | sed 's/^/skipped tests: /'
echo "the rest of the suite (tests that need inputs this checkout lacks fail on a missing file):"
$R pytest -q -p no:cacheprovider -rfE tests --ignore tests/test_rna_genes.py --ignore tests/test_rna_rates.py --ignore tests/test_rna_tests.py > slurm/logs/$SLURM_JOB_ID.rest.log 2>&1
tail -1 slurm/logs/$SLURM_JOB_ID.rest.log
grep -E "^(FAILED|ERROR)" slurm/logs/$SLURM_JOB_ID.rest.log | sed 's/::.*//; s/ - .*//' | sort | uniq -c
grep -E "^E +[A-Za-z]*(Error|Exit)" slurm/logs/$SLURM_JOB_ID.rest.log | sed 's/(os error 2)//; s/[^ ]*\/data\//data\//' | cut -c1-140 | sort | uniq -c | sort -rn | head -12
$R python slurm/h37-mutations.py; b=$?
echo "breaks: exit $b"
[ $t -eq 0 ] && [ $b -eq 0 ]
