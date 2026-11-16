#!/bin/bash
# H17, CPU part: the co-evolution score of the published cybrids (`mitoport cybrid-coev`), with H16's baseline again
# (`mitoport cybrid`), the tests, the breaks of each guard and a rerun check. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision 0f4e970554b2fe0265e718321f422347293dde17
#   uv run --frozen mitoage download pdb 5XTD 5XTE 5Z62 8H9S
#   uv run --frozen mitoage download cybrid
#   the contact maps and the vertebrate partner alignments of H14's job (slurm/job1-contacts.sbatch), copied on the cluster:
#     J=/data/uq2347/runs/mitoage/job1-contacts/data/interim; mkdir -p data/interim/mitoport/partners data/interim/tree
#     cp -r $J/mitoport/iface $J/mitoport/iface-I $J/mitoport/iface-III $J/mitoport/iface-V data/interim/mitoport/
#     cp -r $J/mitoport/partners/aln data/interim/mitoport/partners/
#   mkdir -p slurm/logs && sbatch slurm/h17-cybrid.sbatch
#SBATCH --job-name=mitoage-h17-cybrid
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
sha256sum data/input/mitoage/species.parquet data/input/cybrid/records.gb data/input/pdb/{5XTD,5XTE,5Z62,8H9S}.cif
sha256sum data/interim/mitoport/iface*/{contacts,sites,control_pairs,chains}.tsv
cat data/interim/mitoport/partners/aln/*.afa | sha256sum
O=data/output/mitoport/cybrid
$M mitoport cybrid
$M mitoport cybrid-coev
cp $O/coev_scores.tsv $O/coev_scores.first.tsv
$M mitoport cybrid-coev > /dev/null
if cmp -s $O/coev_scores.tsv $O/coev_scores.first.tsv; then echo "rerun: identical bytes"; else echo "rerun: DIFFERS"; fi
rm $O/coev_scores.first.tsv
$M mitoport cybrid-readout
sha256sum $O/*.tsv
$UV run --frozen pytest -q -p no:cacheprovider tests/test_cybrid.py tests/test_cybrid_h17.py
$UV run --frozen python slurm/h17-cybrid-mutations.py
