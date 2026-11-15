#!/bin/bash
# Every analysis whose numbers the MitoAge manuscript reports (H1 aa-pgls in four classes, H2's column screen in
# mammals and birds, H7 rates, H8 motif-pgls, H11 generation-time in both maturity modes, H12 structure-pgls), on
# one build, in one job, on the cluster host haumea.
# The cluster runs only code that is on main. On the login node, from the root of a checkout pulled from main,
# with MITOAGE_PIPELINES_DATA pointing at a data root of its own for each build (network is needed there, not
# in the job):
#   uv sync --frozen
#   export MITOAGE_PIPELINES_DATA=<root>
#   uv run mitoage download build --revision <dataset commit>
#     (or copy public/data of a dna-seq/mitoage commit into <root>/input/mitoage, for a build not yet on Hugging Face)
#   uv run mitoage download tree --ref <the dna-seq/mitoage commit of that build>
#   uv run mitoage download pdb 5Z62
#   mkdir -p slurm/logs && sbatch --export=ALL,MITOAGE_PIPELINES_DATA=<root> slurm/manuscript-rerun.sbatch
# Writes the tables under <root>/output/longevity and slurm/logs/<job>.log, which records the commit, the inputs'
# digests and the tables' digests. `slurm/compare-runs.py <root A> <root B>` then says, table by table and class
# by class, which rows two such runs share.
#SBATCH --job-name=mitoage-rerun
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
: "${MITOAGE_PIPELINES_DATA:?set MITOAGE_PIPELINES_DATA to the data root of the build to run on}"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
R="$MITOAGE_PIPELINES_DATA"
# git is not installed on the compute nodes: read the commit from .git directly.
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "data root $R"
sha256sum "$R/input/mitoage/species.parquet" "$R/input/mitoage/mtdna_structure.parquet" \
  "$R/input/mitoage/protein_gene_features.parquet" "$R/input/mitoage-site/phylogeny.json" "$R/input/pdb/5Z62.cif"
step() { echo "== $* ($(date -u +%FT%TZ))"; $M "$@"; }
step longevity dated-tree
step longevity controls > /dev/null
step longevity fourfold
for c in Mammalia Aves Actinopterygii Reptilia; do step longevity aa-pgls --class "$c"; done
for c in Mammalia Aves; do step longevity pgls --class "$c"; done
step longevity rates
step longevity motif-pgls
step mitoport contacts
step mitoport burial
step longevity structure-pgls
for m in female any; do step longevity generation-time --maturity "$m"; done
echo "== done ($(date -u +%FT%TZ))"
sha256sum "$R"/output/longevity/*.tsv
