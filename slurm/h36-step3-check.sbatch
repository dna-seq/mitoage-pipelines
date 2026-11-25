#!/bin/bash
# H36, step 3, after the two runs of slurm/h36-step3.sbatch (each host, each twice): collect each host's tasks into the tables and apply
# the step-2 rule, compare every task table and every collected table of the second run with the first byte for byte, the bar's
# comparison with step 1 and step 2, the tests and the breaks of each guard. CPU only (a login-node-sized job on the compute nodes).
#   sbatch --account=<account> --nice=100 --comment '...' --dependency=afterany:<jobs> \
#     --export=ALL,H36_HUMAN=<folder>,H36_MOUSE=<folder>,H36_HUMAN2=<folder>,H36_MOUSE2=<folder> slurm/h36-step3-check.sbatch
#SBATCH --job-name=mitoage-h36-step3-check
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-saprot
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
R="$UV run --frozen --no-sync"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
$R mitoage mitoport saprot-collect --host homo-sapiens --out "$H36_HUMAN"
$R mitoage mitoport saprot-collect --host mus-musculus --out "$H36_MOUSE"
$R mitoage mitoport saprot-compare --human "$H36_HUMAN" --mouse "$H36_MOUSE" --out "$H36_HUMAN"
sha256sum "$H36_HUMAN"/saprot_*.tsv "$H36_MOUSE"/saprot_*.tsv
if [ -n "$H36_HUMAN2" ]; then
  $R mitoage mitoport saprot-collect --host homo-sapiens --out "$H36_HUMAN2" > /dev/null
  $R mitoage mitoport saprot-collect --host mus-musculus --out "$H36_MOUSE2" > /dev/null
  $R mitoage mitoport saprot-compare --human "$H36_HUMAN2" --mouse "$H36_MOUSE2" --out "$H36_HUMAN2" > /dev/null
  for pair in "$H36_HUMAN:$H36_HUMAN2" "$H36_MOUSE:$H36_MOUSE2"; do
    a=${pair%%:*}; b=${pair##*:}
    same=0; differ=0; absent=0
    for f in "$a"/saprot/*.tsv; do
      g="$b/saprot/$(basename $f)"
      if [ ! -f "$g" ]; then absent=$((absent+1)); elif cmp -s "$f" "$g"; then same=$((same+1)); else differ=$((differ+1)); echo "differs: $(basename $f)"; fi
    done
    echo "$a: $same task tables identical byte for byte, $differ differ, $absent not run a second time"
    for f in "$a"/saprot_*.tsv; do
      if cmp -s "$f" "$b/$(basename $f)"; then echo "identical: $(basename $f)"; else echo "differs: $(basename $f)"; fi
    done
  done
else
  echo "second run: not compared (no folder given)"
fi
$R pytest -q -p no:cacheprovider tests/test_cybrid_saprot.py tests/test_cybrid_swap.py
$R python slurm/h36-step3-mutations.py
