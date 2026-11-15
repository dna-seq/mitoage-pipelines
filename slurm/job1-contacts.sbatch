#!/bin/bash
# H14: the contact maps of complexes I (5XTD), III (5XTE) and V (8H9S), with complex IV (5Z62) again, and the readouts
# fixed in advance (`mitoport readouts`), on haumea, with the tests, the mutation checks, a rerun check and a check that
# complex IV's tables equal those of the code before the generalisation. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision 0f4e970554b2fe0265e718321f422347293dde17
#   uv run --frozen mitoage download tree --ref 03c09498564c17d1445bc4c18c70d0e710d471db     # the build's tree (site 03c0949)
#   uv run --frozen mitoage download pdb 5XTD 5Z62 8H9S 5XTE
#   uv run --frozen mitoage download orthologs                                             # H5's 62 genes, mammals
#   uv run --frozen mitoage download orthologs $(printf -- '--gene %s ' $(uv run --frozen python -c \
#     'from mitoage_pipelines.genes import *; print(*COMPLEX_I_NUCLEAR, *COMPLEX_V_NUCLEAR)'))      # 52 more, mammals
#   uv run --frozen mitoage download orthologs --ortholog-taxon vertebrata --dest data/input/ncbi-vertebrata \
#     $(printf -- '--gene %s ' $VERTEBRATA)                    # the subunits that touch, listed in VERTEBRATA below
#   git archive <the commit before the generalisation> | tar -x -C ../job1-contacts-base     # the code compared with
#   mkdir -p slurm/logs && sbatch slurm/job1-contacts.sbatch
# MAFFT 7.526 is taken from $MAFFT (on haumea: the mitoage-msa environment).
#SBATCH --job-name=mitoage-job1-contacts
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
BASE_SRC=${BASE_SRC:-$(cd .. && pwd)/job1-contacts-base/src}
base() {        # the earlier code, reading and writing this checkout's data folder
  MITOAGE_PIPELINES_DATA="$PWD/data" PYTHONPATH="$BASE_SRC" .venv/bin/python -c \
    'import sys; from mitoage_pipelines.cli import main; sys.argv[0] = "mitoage"; main()' "$@"
}
export MAFFT=${MAFFT:-/data/uq2347/micromamba/envs/mitoage-msa/bin/mafft}
VERTEBRATA="COX4I1 COX5A COX5B COX6A1 COX6B1 COX6C COX7A2 COX7B COX7C COX8A NDUFA4 CYC1 UQCR11 UQCRB UQCRC1 UQCRFS1 UQCRQ
NDUFA1 NDUFA3 NDUFA6 NDUFA8 NDUFA9 NDUFA10 NDUFA11 NDUFA12 NDUFA13 NDUFAB1 NDUFB1 NDUFB2 NDUFB3 NDUFB4 NDUFB5 NDUFB6
NDUFB7 NDUFB8 NDUFB9 NDUFB10 NDUFB11 NDUFC2 NDUFS2 NDUFS3 NDUFS5 NDUFS7 NDUFS8 ATP5MC1 ATP5MF ATP5MJ ATP5PB ATP5PD ATP5PF"
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); compared with the code in $BASE_SRC"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet; MAFFT: $($MAFFT --version 2>&1 | head -1)"
D=data; W=data/interim/job1-check; rm -rf $W; mkdir -p $W
$M longevity dated-tree
$M longevity controls > /dev/null
$M longevity fourfold
for c in IV I III V; do
  $M mitoport contacts --complex $c
  $M mitoport burial --complex $c
done
$M mitoport select                                           # mammals, H5's design (data/interim/mitoport/sel)
$M mitoport iface-align --threads 1 --jobs 16                # data/interim/mitoport/iface_aln
P=data/interim/mitoport/partners
$M mitoport select --ncbi data/input/ncbi-vertebrata --out $P/sel --class Mammalia --class Aves --class Actinopterygii --class Reptilia
$M mitoport iface-align --sel $P/sel --out $P/aln --threads 1 --jobs 16
sha256sum data/interim/tree/dated.nwk data/output/longevity/composition.tsv data/output/longevity/fourfold.tsv \
  data/input/pdb/*.cif data/interim/mitoport/iface*/*.tsv data/interim/mitoport/iface_aln/*.afa $P/aln/*.afa

# complex IV through the code before the generalisation, on the same inputs: the tables must be the same bytes
mkdir -p $W/aln62 $W/base-map
for g in $($UV run --frozen python -c 'from mitoage_pipelines.genes import NUCLEAR; print(*NUCLEAR)'); do
  if [ -f data/interim/mitoport/iface_aln/$g.afa ]; then ln -s "$PWD/data/interim/mitoport/iface_aln/$g.afa" $W/aln62/; fi
done
CIF=$D/input/pdb/5Z62.cif
base mitoport contacts --cif $CIF --build $D/input/mitoage --work $W/base-map > /dev/null
base mitoport burial --cif $CIF --build $D/input/mitoage --work $W/base-map > /dev/null
base mitoport iface --skip-esm --work $W/base-map --nucaln $W/aln62 --cif $CIF --build $D/input/mitoage --out $W/base-iface > /dev/null
base mitoport partners --cif $CIF --build $D/input/mitoage --work $W/base-map --nucaln $P/aln --out $W/base-partners > /dev/null
$M mitoport iface --skip-esm --nucaln $W/aln62 --out $W/new-iface62 > /dev/null
same() { if cmp -s "$1" "$2"; then echo "complex IV $3: identical bytes"; else echo "complex IV $3: DIFFERS"; FAIL=1; fi; }
FAIL=0
for f in contacts.tsv sites.tsv control_pairs.tsv chains.tsv species.txt; do
  same $W/base-map/$f data/interim/mitoport/iface/$f "$f, before and after the generalisation"; done
for f in coupling.tsv lifespan.tsv; do same $W/base-iface/$f $W/new-iface62/$f "$f on H5's 62 alignments, before and after"; done

$UV run --frozen pytest -q -p no:cacheprovider tests/test_complexes.py tests/test_partners.py tests/test_iface_structure.py
$UV run --frozen python slurm/job1-contacts-mutations.py
$UV run --frozen python slurm/structure-pgls-mutations.py cpu

steps() {
  for c in IV I III V; do
    $M mitoport iface --complex $c --skip-esm
    $M mitoport partners --complex $c
  done
  $M mitoport readouts
}
OUTS="data/output/mitoport/iface/coupling.tsv data/output/mitoport/iface/lifespan.tsv data/output/mitoport/partners.tsv
data/output/mitoport/partners_species.tsv data/output/mitoport/complex-*/*.tsv data/output/mitoport/contact_readouts.tsv"
steps
for f in coupling.tsv lifespan.tsv; do same $W/new-iface62/$f data/output/mitoport/iface/$f "$f on every alignment against H5's 62 only"; done
for f in partners.tsv partners_species.tsv; do same $W/base-partners/$f data/output/mitoport/$f "$f, before and after"; done
mkdir -p $W/first; for f in $OUTS; do mkdir -p $W/first/$(dirname $f); cp $f $W/first/$f; done
steps > /dev/null
for f in $OUTS; do if cmp -s $f $W/first/$f; then echo "rerun $f: identical bytes"; else echo "rerun $f: DIFFERS"; FAIL=1; fi; done
sha256sum $OUTS
echo "checks failed: $FAIL"
exit $FAIL
