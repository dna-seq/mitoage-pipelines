#!/bin/bash
# H38: the tests of the architecture code and the breaks of its guards (slurm/h38-mutations.py), then the rest of the suite, from a fresh
# clone of main after slurm/h38-run.sbatch (the tests read the tables it made). H38_TOOLS names the folder that holds Tandem Repeats
# Finder. CPU only.
#   mkdir -p slurm/logs && sbatch --account=<account> --partition=compute --nice=100 --export=ALL,H38_TOOLS=<folder> --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h38-check.sbatch
#SBATCH --job-name=mitoage-h38-check
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
R="uv run --frozen"
set +e
$R pytest -q -p no:cacheprovider -rs tests/test_arch.py; t=$?
echo "the tests of H38: exit $t (a skipped test is listed above and counts as a failure of this check)"
$R pytest -q -p no:cacheprovider -rs tests/test_arch.py 2>&1 | grep -c "^SKIPPED" | sed 's/^/skipped tests: /'
echo "the rest of the suite (tests that need inputs this checkout lacks fail on a missing file):"
$R pytest -q -p no:cacheprovider -rfE tests --ignore tests/test_arch.py > slurm/logs/$SLURM_JOB_ID.rest.log 2>&1
tail -1 slurm/logs/$SLURM_JOB_ID.rest.log
grep -E "^(FAILED|ERROR)" slurm/logs/$SLURM_JOB_ID.rest.log | sed 's/::.*//; s/ - .*//' | sort | uniq -c
$R python slurm/h38-mutations.py; b=$?
echo "breaks: exit $b"
[ $t -eq 0 ] && [ $b -eq 0 ]
