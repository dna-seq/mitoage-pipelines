#!/bin/bash
# H26's R1 (HyPhy RELAX, --models Minimal) on the long-lived clades (draw 0) and 20 matched control draws, and R4
# (aBSREL on the stems) on draws 0 and 1, for mammals and birds; one array task per run. After
# slurm/adapt-longevity.sbatch, in the same clone:
#   sbatch --dependency=afterok:<adapt-longevity job> slurm/adapt-longevity-selection.sbatch
# Tasks 0-41: RELAX, class = CLASSES[task / 21], draw = task % 21. Tasks 42-45: aBSREL, class = CLASSES[(task-42) / 2],
# draw = (task-42) % 2. HyPhy 2.5.104: set HYPHY to its binary, or have `hyphy` on PATH. The Slurm account is given on
# the command line (`sbatch --account=...`). No time limit is set: no RELAX run of this size has been timed yet; the
# first task's log gives the rate.
#SBATCH --job-name=mitoage-adapt-longevity-selection
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-45
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
if [ -n "$HYPHY" ]; then export PATH="$(dirname "$HYPHY"):$PATH"; fi
command -v hyphy >/dev/null || { echo "hyphy is not on PATH: set HYPHY to its binary"; exit 1; }
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
CLASSES=(Mammalia Aves)
O=data/output/adapt-longevity
T=$SLURM_ARRAY_TASK_ID
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); $(hyphy --version 2>&1 | head -1)"
if [ "$T" -lt 42 ]; then
  cls=${CLASSES[$((T / 21))]}; draw=$((T % 21))
  echo "RELAX long_lived $cls draw $draw"
  $M adapt relax --trait long_lived --draw $draw --cpus 4 --out $O/$cls
else
  k=$((T - 42)); cls=${CLASSES[$((k / 2))]}; draw=$((k % 2))
  echo "aBSREL long_lived $cls draw $draw"
  $M adapt absrel --trait long_lived --draw $draw --cpus 4 --out $O/$cls
fi
