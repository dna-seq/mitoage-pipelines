#!/bin/bash
# H21-H23 R1 and H25 B4 (HyPhy RELAX, --models Minimal) on the trait clades (draw 0) and 20 matched control draws, and
# R4 (aBSREL on the stems) on draws 0 and 1; one array task per run. After slurm/adapt-main.sbatch, in the same clone:
#   sbatch --dependency=afterok:<adapt-main job> slurm/adapt-selection.sbatch
# Tasks 0-104: RELAX, trait = TRAITS[task / 21], draw = task % 21. Tasks 105-114: aBSREL, trait = TRAITS[(task-105) / 2],
# draw = (task-105) % 2. HyPhy 2.5.104 from $HYPHY (on haumea: the mitoage-msa environment).
#SBATCH --job-name=mitoage-adapt-selection
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --array=0-114
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
export PATH="$(dirname ${HYPHY:-/data/uq2347/micromamba/envs/mitoage-msa/bin/hyphy}):$PATH"
TRAITS=(diadromy deep_habitat marine hibernation bat_hibernation)
T=$SLURM_ARRAY_TASK_ID
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); $(hyphy --version 2>&1 | head -1)"
if [ "$T" -lt 105 ]; then
  trait=${TRAITS[$((T / 21))]}; draw=$((T % 21))
  echo "RELAX $trait draw $draw"
  $M adapt relax --trait $trait --draw $draw --cpus 4
else
  k=$((T - 105)); trait=${TRAITS[$((k / 2))]}; draw=$((k % 2))
  echo "aBSREL $trait draw $draw"
  $M adapt absrel --trait $trait --draw $draw --cpus 4
fi
