#!/bin/bash
# H15 (RM16): the published mtDNA-encoded peptides across the build, the window, synonymous and dinucleotide nulls, the
# tree and lifespan (`mitoage peptides ...`), on haumea, with the tests, the mutation check and a rerun of every step.
# The cluster runs only code that is on main. On the login node, from the root of a checkout pulled from main
# (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision c525adbff8c12b19981fdada7f1615c2559b7140
#   uv run --frozen mitoage download tree --ref 6fe8e726fb32dd1dddc4a967db72af4f60c0215c
#   mkdir -p slurm/logs && sbatch slurm/job4-peptides.sbatch
# The tree ref is the one H12's runs downloaded; the site's newer trees name species build c525adb does not have.
#SBATCH --job-name=mitoage-job4-peptides
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
# git is not installed on the compute nodes: read the commit from .git directly.
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet; tree: $(sha256sum data/input/mitoage-site/phylogeny.json | cut -c1-12) phylogeny.json"
$M longevity dated-tree
$M longevity controls > /dev/null
$M longevity fourfold
sha256sum data/interim/tree/dated.nwk data/output/longevity/composition.tsv data/output/longevity/fourfold.tsv
O=data/output/peptides
run_all() {
  $M peptides catalogue
  $M peptides frames
  $M peptides nulls --jobs 8
  $M peptides tree
  $M peptides pgls
}
time run_all
TABLES="$O/catalogue.tsv $O/left_out.tsv $O/states.tsv $O/kept.tsv $O/nulls.tsv $O/tree.tsv $O/pgls.tsv"
sha256sum $TABLES
$UV run --frozen pytest -q -p no:cacheprovider tests/test_peptides.py
$UV run --frozen python slurm/job4-peptides-mutations.py
mkdir -p $O/first && cp $TABLES $O/first/
run_all > /dev/null
for t in $TABLES; do
  if cmp -s "$t" "$O/first/$(basename $t)"; then echo "rerun $(basename $t): identical bytes"; else echo "rerun $(basename $t): DIFFERS"; fi
done
rm -r $O/first
sha256sum $TABLES
