#!/bin/bash
# H19: mito-nuclear coupling in mammals, birds and ray-finned fishes, its fall with distance from the mtDNA-encoded
# core, and its mitochondrial side (`mitoport classes`), on haumea, with the tests, the breaks of each guard and a
# rerun check. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision 0f4e970554b2fe0265e718321f422347293dde17
#   uv run --frozen mitoage download tree --ref 03c09498564c17d1445bc4c18c70d0e710d471db     # the build's tree (site 03c0949)
#   uv run --frozen mitoage download pdb 5XTD 5Z62 8H9S 5XTE
#   # the 50 vertebrate ortholog sets of H14's job, as it downloaded them (NCBI's sets change with every annotation
#   # release; their records are in that clone's data/input/MANIFEST.json):
#   cp -r /data/uq2347/runs/mitoage/job1-contacts/data/input/ncbi-vertebrata data/input/
#   # the other genes of H19 (a set already there is skipped, and said so):
#   uv run --frozen mitoage download orthologs --ortholog-taxon vertebrata --dest data/input/ncbi-vertebrata \
#     $(printf -- '--gene %s ' $(uv run --frozen python -c 'from mitoage_pipelines.genes import H19_NUCLEAR; print(*H19_NUCLEAR)'))
#   mkdir -p slurm/logs && sbatch slurm/h19-classes.sbatch
# MAFFT 7.526 is taken from $MAFFT (on haumea: the mitoage-msa environment).
#SBATCH --job-name=mitoage-h19-classes
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
export MAFFT=${MAFFT:-/data/uq2347/micromamba/envs/mitoage-msa/bin/mafft}
H14_ALN=${H14_ALN:-/data/uq2347/runs/mitoage/job1-contacts/data/interim/mitoport/partners/aln}
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet; MAFFT: $($MAFFT --version 2>&1 | head -1)"
echo "vertebrate ortholog sets: $(ls data/input/ncbi-vertebrata | wc -l)"
FAIL=0
$M longevity dated-tree
# the tests and the breaks of each guard first: they need the build, the tree and the structures only
$UV run --frozen pytest -q -p no:cacheprovider tests/test_classes.py
$UV run --frozen python slurm/h19-classes-mutations.py
for c in IV I III V; do
  $M mitoport contacts --complex $c
  $M mitoport burial --complex $c
done
P=data/interim/mitoport/partners
$M mitoport select --ncbi data/input/ncbi-vertebrata --out $P/sel --class Mammalia --class Aves --class Actinopterygii --class Reptilia
$M mitoport iface-align --sel $P/sel --out $P/aln --threads 1 --jobs 16
# the 50 alignments H14 read, made again here from the same ortholog sets: the same bytes?
if [ -d "$H14_ALN" ]; then
  same=0; differ=0
  for f in "$H14_ALN"/*.afa; do
    if cmp -s "$f" $P/aln/$(basename "$f"); then same=$((same + 1)); else differ=$((differ + 1)); echo "differs from H14's alignment: $(basename "$f")"; fi
  done
  echo "alignments of H14's job made again: $same identical, $differ different"
else
  echo "H14's alignments are not at $H14_ALN: NOT COMPARED"
fi
sha256sum data/interim/tree/dated.nwk data/input/pdb/*.cif data/interim/mitoport/iface*/{contacts,sites,control_pairs,chains}.tsv
cat $P/aln/*.afa | sha256sum
O=data/output/mitoport/classes
$M mitoport classes
mkdir -p $O.first && cp $O/*.tsv $O.first/
$M mitoport classes > /dev/null
for f in $O/*.tsv; do
  if cmp -s $f $O.first/$(basename $f); then echo "rerun $(basename $f): identical bytes"; else echo "rerun $(basename $f): DIFFERS"; FAIL=1; fi
done
rm -r $O.first
sha256sum $O/*.tsv
echo "checks failed: $FAIL"
exit $FAIL
