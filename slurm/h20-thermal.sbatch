#!/bin/bash
# H20: temperature, size and metabolic rate against the mtDNA substitution rate of sister pairs and its composition,
# and Lehmann et al.'s GC and metabolic rate retested on the tree (`longevity thermal`), on haumea, with the tests,
# the breaks of each guard and a rerun check. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision 0f4e970554b2fe0265e718321f422347293dde17
#   uv run --frozen mitoage download tree --ref 03c09498564c17d1445bc4c18c70d0e710d471db     # the build's tree (site 03c0949)
#   uv run --frozen mitoage download fishbase-temperature                                    # md5-checked
#   mkdir -p slurm/logs && sbatch slurm/h20-thermal.sbatch
#SBATCH --job-name=mitoage-h20-thermal
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet"
FB=data/input/build2/fishbase/v26.06
md5sum $FB/species.parquet $FB/synonyms.parquet $FB/stocks.parquet
FAIL=0
$M longevity dated-tree
$M longevity controls > /dev/null
$M longevity fourfold
sha256sum data/interim/tree/dated.nwk data/output/longevity/fourfold.tsv
$UV run --frozen pytest -q -p no:cacheprovider tests/test_thermal.py
$UV run --frozen python slurm/h20-thermal-mutations.py
O=data/output/longevity/thermal
$M longevity thermal
mkdir -p $O.first && cp $O/*.tsv $O.first/
$M longevity thermal > /dev/null
for f in $O/*.tsv; do
  if cmp -s $f $O.first/$(basename $f); then echo "rerun $(basename $f): identical bytes"; else echo "rerun $(basename $f): DIFFERS"; FAIL=1; fi
done
rm -r $O.first
sha256sum $O/*.tsv
echo "checks failed: $FAIL"
exit $FAIL
