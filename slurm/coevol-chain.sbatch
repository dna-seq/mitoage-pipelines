#!/bin/bash
# One 3-day leg of a Coevol chain for H3 (`longevity coevol-resume`), on the cluster host haumea. Submit legs of
# one chain with --dependency=afterany on the previous leg; each continues from the chain's last saved point
# (checkpoint every 6 h) and stops starting segments after 70 h, which leaves the 72 h limit room for one
# cycle (about 65-80 s), the 10-minute guard that ends Coevol's MPI build, and the checkpoint.
#   CHAIN=1 OUT=<data root>/output/longevity/coevol            the published chains, copied there first
#   CHAIN=1 OUT=<data root>/output/longevity/coevol-fixed FIXED=1   chains holding the dated ages (-fixtimes);
#                                                 the first leg starts the chain from coevol-prep --fixed-ages
#   MPI=8   ranks of Coevol's MPI build (Makefile.mpi of Coevol 1.6, built with gcc 15.2.0 and OpenMPI 5.0.8,
#           at $COEVOL_MPI); MPI=0 runs the serial build at $COEVOL.
# The cluster runs only code that is on main: from the root of a checkout pulled from main, with
# MITOAGE_PIPELINES_DATA set,  mkdir -p slurm/logs && sbatch --export=ALL,CHAIN=1,OUT=...,MPI=8 slurm/coevol-chain.sbatch
#SBATCH --job-name=mitoage-coevol
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm/logs/coevol-%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${MITOAGE_PIPELINES_DATA:?}" "${CHAIN:?}" "${OUT:?}"
MPI=${MPI:-8}; UNTIL=${UNTIL:-15000}
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); job $SLURM_JOB_ID on $(hostname); chain $CHAIN in $OUT; MPI $MPI; until $UNTIL"
module load gcc/15.2.0 openmpi/gcc.15/5.0.8
if [ "$MPI" -gt 0 ]; then export COEVOL=${COEVOL_MPI:-/data/uq2347/tools/coevol-mpi/data/coevol}; else export COEVOL=${COEVOL:-/data/uq2347/tools/coevol/data/coevol}; fi
M="${UV:-uv} run --frozen mitoage"
if [ "${FIXED:-0}" = 1 ] && [ ! -e "$OUT/chain$CHAIN.param" ]; then
  # the first point with the serial build (one cycle), then the legs continue it
  COEVOL=/data/uq2347/tools/coevol/data/coevol $M longevity coevol --fixed-ages --chain "$CHAIN" --cycles 1 --out "$OUT"
fi
$M longevity coevol-resume --chain "$CHAIN" --until "$UNTIL" --stop-after-hours 70 --mpi "$MPI" --out "$OUT"
echo "== leg done $(date -u +%FT%TZ)"; tail -1 "$OUT/chain$CHAIN.trace" | cut -f1-4
