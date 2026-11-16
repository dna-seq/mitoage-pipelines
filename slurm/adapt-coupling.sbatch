#!/bin/bash
# H21-H23 C (and H23's positive control P, the bats' stem): per-protein branch lengths on the dated topology pruned to the
# species with every protein, then C1, C2, P. After slurm/adapt-main.sbatch, in the same clone. The nuclear alignments
# are H19's (vertebrate ortholog sets, one protein per species, MAFFT L-INS-i on one thread): pass H19_ALN=<folder> when
# H19's job has made them; otherwise they are made here by H19's own commands, which need the ortholog sets on the login
# node first, as in slurm/h19-classes.sbatch:
#   uv run --frozen mitoage download orthologs --ortholog-taxon vertebrata --dest data/input/ncbi-vertebrata \
#     $(printf -- '--gene %s ' $(uv run --frozen python -c 'from mitoage_pipelines.genes import H19_NUCLEAR; print(*H19_NUCLEAR)'))
#   sbatch --dependency=afterok:<adapt-main job> slurm/adapt-coupling.sbatch
# Then, when this and slurm/adapt-selection.sbatch are done: uv run --frozen mitoage adapt readout
# IQ-TREE 3.1.4 and MAFFT 7.526 from the mitoage-msa environment.
#SBATCH --job-name=mitoage-adapt-coupling
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
ENV=${ENV:-/data/uq2347/micromamba/envs/mitoage-msa/bin}
export PATH="$ENV:$PATH" MAFFT=${MAFFT:-$ENV/mafft}
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); $(iqtree --version 2>&1 | head -1)"
P=data/interim/mitoport/partners
if [ -n "$H19_ALN" ]; then
  ALN=$H19_ALN; echo "nuclear alignments: H19's, $ALN"
else
  $M mitoport select --ncbi data/input/ncbi-vertebrata --out $P/sel --class Mammalia --class Aves --class Actinopterygii --class Reptilia
  $M mitoport iface-align --sel $P/sel --out $P/aln --threads 1 --jobs 16
  ALN=$P/aln
fi
cat $ALN/*.afa | sha256sum
O=data/output/adaptations
for cls in Mammalia Actinopterygii; do
  W=$O/coupling/$cls
  $M adapt coupling-prep --class $cls --nucaln $ALN
  cat $W/prepare.tsv
  $M mitoport branch-lengths --aln-dir $W --topology $W/topology.nwk --genes $W/genes.txt --bl-dir $W/bl --threads 16 --seed 1
  $M adapt coupling --class $cls
done
