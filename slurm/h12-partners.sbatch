#!/bin/bash
# H12, second round, CPU part: hydropathy with each gene's own composition (`longevity structure-pgls`) and the
# discordance of complex IV's mito-nuclear contact pairs (`mitoport partners`), on haumea, with the tests, the
# mutation check and rerun checks. The cluster runs only code that is on main. On the login node, from the root of a
# checkout pulled from main (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision c525adbff8c12b19981fdada7f1615c2559b7140
#   uv run --frozen mitoage download tree && uv run --frozen mitoage download pdb 5Z62
#   uv run --frozen mitoage download orthologs --ortholog-taxon vertebrata --dest data/input/ncbi-vertebrata \
#     --gene COX4I1 --gene COX5A --gene COX5B --gene COX6A1 --gene COX6B1 --gene COX6C --gene COX7A2 --gene COX7B \
#     --gene COX7C --gene COX8A --gene NDUFA4
#   mkdir -p slurm/logs && sbatch slurm/h12-partners.sbatch
# MAFFT 7.526 is taken from $MAFFT (on haumea: the mitoage-msa environment).
#SBATCH --job-name=mitoage-h12-partners
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
export MAFFT=${MAFFT:-/data/uq2347/micromamba/envs/mitoage-msa/bin/mafft}
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet; MAFFT: $($MAFFT --version 2>&1 | head -1)"
$M longevity dated-tree
$M longevity controls > /dev/null
$M longevity fourfold
$M mitoport contacts
$M mitoport burial
P=data/interim/mitoport/partners
$M mitoport select --ncbi data/input/ncbi-vertebrata --out $P/sel --class Mammalia --class Aves --class Actinopterygii --class Reptilia
$M mitoport iface-align --sel $P/sel --out $P/aln --threads 1 --jobs 8
sha256sum data/interim/tree/dated.nwk data/output/longevity/composition.tsv data/output/longevity/fourfold.tsv \
  data/interim/mitoport/iface/sites.tsv $P/aln/*.afa
$UV run --frozen pytest -q -p no:cacheprovider tests/test_structure_pgls.py tests/test_partners.py tests/test_motif_pgls.py
$UV run --frozen python slurm/structure-pgls-mutations.py cpu
for step in "longevity structure-pgls:data/output/longevity/structure_pgls.tsv" "mitoport partners:data/output/mitoport/partners.tsv"; do
  cmd=${step%%:*}; table=${step#*:}
  $M $cmd
  cp "$table" "$table.first"
  $M $cmd > /dev/null
  if cmp -s "$table" "$table.first"; then echo "rerun $cmd: identical bytes"; else echo "rerun $cmd: DIFFERS"; fi
  rm "$table.first"
done
sha256sum data/output/longevity/structure_pgls.tsv data/output/mitoport/partners.tsv data/output/mitoport/partners_species.tsv
