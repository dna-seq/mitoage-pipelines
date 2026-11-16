#!/bin/bash
# H26 (docs/HYPOTHESES.md), first job: the tests of the trait, the breaks of each guard, the trait states, the pairs
# per origin and the RELAX/aBSREL inputs for mammals and birds; then the same steps again into a second folder, to
# check that they come out the same bytes. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision 0f4e970554b2fe0265e718321f422347293dde17
#   uv run --frozen mitoage download tree --ref 03c09498564c17d1445bc4c18c70d0e710d471db     # the build's tree (site 03c0949)
#   mkdir -p slurm/logs && sbatch --account=<the Slurm account> slurm/adapt-longevity.sbatch
# Then slurm/adapt-longevity-selection.sbatch (RELAX, aBSREL); `adapt longevity-coupling --work <H23's coupling/Mammalia>`
# once H23's branch lengths exist; `adapt longevity-readout` last.
#SBATCH --job-name=mitoage-adapt-longevity
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 POLARS_MAX_THREADS=4
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet"
$M longevity dated-tree
$UV run --frozen pytest -q -p no:cacheprovider tests/test_adapt_longevity.py
$UV run --frozen python slurm/adapt-longevity-mutations.py
O=data/output/adapt-longevity
run_all () {
  local out=$1
  $M adapt longevity-traits --out $out
  $M adapt longevity-pairs --out $out
  for c in Mammalia Aves; do
    $M adapt longevity-selection-prep --class $c --out $out
  done
}
run_all $O
run_all $O.rerun
if diff -rq $O $O.rerun; then echo "same bytes: every table and selection input"; FAIL=0; else echo "DIFFERS on rerun"; FAIL=1; fi
rm -rf $O.rerun
sha256sum $O/*.tsv $O/*/*.tsv
exit $FAIL
