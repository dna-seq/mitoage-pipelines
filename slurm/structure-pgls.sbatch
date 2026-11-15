#!/bin/bash
# H12 (`mitoage longevity structure-pgls`) and H8 (`longevity motif-pgls`) on the cluster host haumea, with the
# tests, a check that each test fails when the code it guards is broken, and a rerun check.
# The cluster runs only code that is on main. On the login node, from the root of a checkout pulled from main
# (network is needed there, not in the job):
#   uv sync --frozen
#   uv run mitoage download build --revision c525adbff8c12b19981fdada7f1615c2559b7140
#   uv run mitoage download tree
#   uv run mitoage download pdb 5Z62
#   mkdir -p slurm/logs && sbatch slurm/structure-pgls.sbatch
# The job makes the remaining inputs itself (dated tree, composition, fourfold counts, the complex IV contact map),
# then writes data/output/longevity/structure_pgls.tsv and motif_pgls.tsv, and slurm/logs/<job>.log, which records
# the commit, the inputs' digests, the tests, the mutation check and whether a second run gives the same bytes.
#SBATCH --job-name=mitoage-h12
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
# git is not installed on the compute nodes: read the commit from .git directly.
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet"
echo "5Z62: $(sha256sum data/input/pdb/5Z62.cif | cut -c1-12)"
$M longevity dated-tree
$M longevity controls > /dev/null          # composition.tsv (and the earlier tree-free tests, not read here)
$M longevity fourfold
$M mitoport contacts
$M mitoport burial
sha256sum data/interim/tree/dated.nwk data/output/longevity/composition.tsv data/output/longevity/fourfold.tsv \
  data/interim/mitoport/iface/sites.tsv
$UV run --frozen pytest -q -p no:cacheprovider tests/test_structure_pgls.py tests/test_motif_pgls.py
$UV run --frozen python slurm/structure-pgls-mutations.py
out=data/output/longevity
$M longevity structure-pgls
cp "$out/structure_pgls.tsv" "$out/structure_pgls.first.tsv"
$M longevity structure-pgls > /dev/null
if cmp -s "$out/structure_pgls.tsv" "$out/structure_pgls.first.tsv"; then echo "rerun: identical bytes"; else echo "rerun: DIFFERS"; fi
rm "$out/structure_pgls.first.tsv"
$M longevity motif-pgls
sha256sum "$out/structure_pgls.tsv" "$out/motif_pgls.tsv"
