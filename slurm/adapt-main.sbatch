#!/bin/bash
# H21-H25 (docs/HYPOTHESES.md, docs/MITO_LOAD_ADAPTATIONS.md), first job: the tests, the breaks of each guard, the trait
# states, the pairs per origin, the phylogenetic fits and the RELAX/aBSREL inputs; then the same steps again into a
# second folder, to check that they come out the same bytes. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision 0f4e970554b2fe0265e718321f422347293dde17
#   uv run --frozen mitoage download tree --ref 03c09498564c17d1445bc4c18c70d0e710d471db     # the build's tree (site 03c0949)
#   uv run --frozen mitoage download fishbase-temperature
#   uv run --frozen mitoage download taxonomy --archive 2026-11-01          # names.dmp: NCBI synonyms for the hibernation tables
#   # the mammal hibernation table (Dryad serves it only to a browser): copy the maintainer's file, then check it
#   mkdir -p data/input/hibernation && cp <the file saved from Dryad> data/input/hibernation/heldstab_2018_supplementary_material.xlsx
#   uv run --frozen mitoage download hibernation                            # fetches the bat table; refuses either file on a wrong digest
#   mkdir -p slurm/logs && sbatch slurm/adapt-main.sbatch
# Then slurm/adapt-selection.sbatch (RELAX, aBSREL) and slurm/adapt-coupling.sbatch, and `adapt readout` last.
#SBATCH --job-name=mitoage-adapt-main
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
NAMES=data/input/build2/taxdump/2026-11-01/names.dmp
export MITOAGE_TAXDUMP=$(dirname $NAMES)
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet"
sha256sum data/input/hibernation/* data/input/build2/fishbase/v26.06/*.parquet $NAMES
$M longevity dated-tree
$M longevity fourfold
# the tests and the breaks of each guard first
$UV run --frozen pytest -q -p no:cacheprovider tests/test_adapt.py tests/test_hibernation.py
$UV run --frozen python slurm/adapt-mutations.py
O=data/output/adaptations
run_all () {
  local out=$1
  $M adapt traits --names-dmp $NAMES --out $out
  $M adapt pairs --out $out
  $M adapt fits --out $out
  for t in diadromy deep_habitat marine hibernation bat_hibernation; do
    $M adapt selection-prep --trait $t --out $out
  done
}
run_all $O
run_all $O.rerun
FAIL=0
for f in traits.tsv pairs.tsv pair_tests.tsv composition.tsv lifespan.tsv l2.tsv bat_pairs.tsv; do
  if cmp -s $O/$f $O.rerun/$f; then echo "same bytes: $f"; else echo "DIFFERS on rerun: $f"; FAIL=1; fi
done
if diff -rq $O/selection $O.rerun/selection; then echo "same bytes: selection/"; else echo "DIFFERS on rerun: selection/"; FAIL=1; fi
rm -rf $O.rerun
sha256sum $O/*.tsv
exit $FAIL
