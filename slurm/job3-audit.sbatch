#!/bin/bash
# H13: the audit of the build's mitogenome records (strand, gene order, genetic code, duplicates and mislabels), on
# haumea, with the tests, the mutation check and a rerun-identical check. The cluster runs only code that is on main.
# On the login node, from the root of a fresh clone of main (network is needed there, not in the job):
#   uv sync --frozen
#   uv run --frozen mitoage download build --revision c525adbff8c12b19981fdada7f1615c2559b7140
#   uv run --frozen mitoage download refseq --from-snapshot     # RefSeq release 237, md5-checked
#   uv run --frozen mitoage download contig --from-snapshot     # the CONTIG-only records' sequence
#   uv run --frozen mitoage download taxonomy --from-snapshot   # NCBI Taxonomy of 2026-11-13 (genetic codes)
#   mkdir -p slurm/logs && sbatch slurm/job3-audit.sbatch
#SBATCH --job-name=mitoage-job3-audit
#SBATCH --partition=compute
#SBATCH --account=itmz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/%j.log
set -eo pipefail
cd "$SLURM_SUBMIT_DIR"
UV=${UV:-uv}
M="$UV run --frozen mitoage"
OUT=data/output/audit
echo "commit $(cat .git/refs/heads/main 2>/dev/null || cat .git/HEAD)"
echo "build: $(sha256sum data/input/mitoage/species.parquet | cut -c1-12) species.parquet"
steps() {
  $M audit records
  $M audit strand
  $M audit order
  $M audit code
  $M audit duplicates
}
steps
$UV run --frozen pytest -q -p no:cacheprovider tests/test_audit.py
$UV run --frozen python slurm/job3-audit-mutations.py
rm -rf "$OUT.first" && cp -r "$OUT" "$OUT.first"
steps > /dev/null 2>&1
if diff -r "$OUT" "$OUT.first" > /dev/null; then echo "rerun of every audit step: identical bytes"; else echo "rerun: DIFFERS"; diff -rq "$OUT" "$OUT.first" || true; fi
rm -rf "$OUT.first"
sha256sum $OUT/*.tsv
