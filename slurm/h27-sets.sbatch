#!/bin/bash
# H27: the manuscript's fast analyses (H1 aa-pgls in four classes and the column screen in mammals and birds, H7
# rates, H8 motif-pgls, H11 generation-time in both maturity modes, H20 thermal) on one build, once unrestricted and
# once on each species set of `longevity/species_sets.py`, as one array task per set. Each task writes into a data
# root of its own, <base>/<set>, whose input is the shared download; nothing of an earlier run is overwritten.
# The cluster runs only code that is on main. On the login node, from the root of a fresh clone of main (network is
# needed there, not in the job), with <base> a folder of its own for this run:
#   uv sync --frozen
#   export MITOAGE_PIPELINES_DATA=<base>/all
#   uv run --frozen mitoage download build --revision <dataset commit of a build that has data_confidence>
#   uv run --frozen mitoage download tree --ref <the dna-seq/mitoage commit of that build>
#   uv run --frozen mitoage download fishbase-temperature
#   for s in not_unreliable high well_sampled captive wild; do mkdir -p <base>/$s && ln -s ../all/input <base>/$s/input; done
#   mkdir -p slurm/logs && sbatch --account=<account> --array=0-5%2 --export=ALL,H27_BASE=<base> slurm/h27-sets.sbatch
# Task 0 is the unrestricted run on this build, the one every set is read against.
#SBATCH --job-name=mitoage-h27-sets
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${H27_BASE:?set H27_BASE to the folder that holds the data roots of this run}"
SETS=(all not_unreliable high well_sampled captive wild)
S=${SETS[$SLURM_ARRAY_TASK_ID]:?no species set for array task $SLURM_ARRAY_TASK_ID}
export MITOAGE_PIPELINES_DATA="$H27_BASE/$S"
if [ "$S" = all ]; then unset MITOAGE_SPECIES_SET; else export MITOAGE_SPECIES_SET="$S"; fi
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
UV=${UV:-uv}
M="$UV run --frozen mitoage"
R="$MITOAGE_PIPELINES_DATA"
# git is not installed on the compute nodes: read the commit from .git directly.
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "species set: ${MITOAGE_SPECIES_SET:-none (all species of each analysis)}"
echo "data root $R"
sha256sum "$R/input/mitoage/species.parquet" "$R/input/mitoage-site/phylogeny.json"
FAIL=0
# One step failing (a set too small for an analysis) must not hide the others: each is run, and its end is said.
step() {
  echo "== $* ($(date -u +%FT%TZ))"
  if $M "$@" > "${OUT:-/dev/stdout}"; then echo "-- ok: $*"; else echo "-- FAILED ($?): $*"; FAIL=1; fi
}
step longevity dated-tree
OUT=/dev/null step longevity controls        # its table is long; the step's end is still said
step longevity fourfold
for c in Mammalia Aves Actinopterygii Reptilia; do step longevity aa-pgls --class "$c"; done
for c in Mammalia Aves; do step longevity pgls --class "$c"; done
step longevity rates
step longevity motif-pgls
for m in female any; do step longevity generation-time --maturity "$m"; done
step longevity thermal
echo "== done, failures: $FAIL ($(date -u +%FT%TZ))"
find "$R/output" -name '*.tsv' | sort | xargs sha256sum
exit $FAIL
