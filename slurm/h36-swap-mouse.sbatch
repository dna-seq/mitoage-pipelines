#!/bin/bash
# H36, step 2 for the mouse host (the swap, `mitoage mitoport swap --host mus-musculus`): StaB-ddG's binding of the donor's
# mitochondrial proteins to one nuclear subunit of the mouse's own solved complex (8OM1, 7O3H, 7O37), with the mouse subunit,
# with the donor's RefSeq ortholog and with random mammals', one array task per (donor, complex, nuclear gene), on CPU nodes; no
# GPU is used. The cluster runs only code that is on main. On the login node, from the root of a fresh clone of main that holds
# the inputs (the structures, the cybrid records, the build, the vertebrate alignments of the nuclear subunits in
# data/interim/mitoport/partners/aln, the weights in data/input/stabddg and the RefSeq proteins in data/input/donor-orthologs; the
# job needs no network):
#   UV_PROJECT_ENVIRONMENT=.venv-stab uv sync --frozen --group stab --group dev
#   mkdir -p slurm/logs && sbatch --account=<account> --nice=100 --comment '<seat=...; what=...; eta=...; cuda=->' slurm/h36-swap-mouse.sbatch
#   sbatch --account=<account> --nice=100 --comment '...' --dependency=afterany:<array job> [--export=ALL,H36_SECOND=<folder>] slurm/h36-swap-mouse-check.sbatch
# Every task a second time into another folder, for the byte comparison of the check job:
#   sbatch --account=<account> --nice=100 --comment '...' --export=ALL,H36_OUT=<folder> slurm/h36-swap-mouse.sbatch
# The tasks are those of `mitoage mitoport swap --host mus-musculus --list` (171: rat, M. caroli and M. pahari, complexes I, III and IV);
# a task whose gene is not assessed ends at once.
#SBATCH --job-name=mitoage-h36-swap-mouse
#SBATCH --partition=compute
#SBATCH --array=0-170%30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=slurm/logs/%A_%a.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
export UV_PROJECT_ENVIRONMENT=.venv-stab
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK POLARS_MAX_THREADS=$SLURM_CPUS_PER_TASK
OUT=${H36_OUT:-data/output/mitoport/cybrid-swap-mouse}
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD); task $SLURM_ARRAY_TASK_ID on $(hostname); out $OUT"
sha256sum data/input/cybrid/records.gb data/input/pdb/{8OM1,7O3H,7O37}.cif data/input/stabddg/*.pt data/input/mitoage/species.parquet data/input/donor-orthologs/*.faa
s=$(date +%s)
# --no-sync: the tasks start together on one environment, and a `uv run` that syncs it removes and rewrites torch's files
# under the others (four of eight tasks of H34 died of it on 2026-11-23)
$UV run --frozen --no-sync mitoage mitoport swap --host mus-musculus --task "$SLURM_ARRAY_TASK_ID" --out "$OUT"
echo "task $SLURM_ARRAY_TASK_ID: $(( $(date +%s) - s )) s"
