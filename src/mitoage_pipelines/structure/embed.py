"""ESM-2 650M embeddings of COX1-3 at complex IV's contact and control residues, for H12's embedding test.

For every species H12 tests (H1's default dataset in each of the four classes), the ungapped COX1, COX2 and COX3
of the build's alignments are run through ESM-2 650M (facebook/esm2_t33_650M_UR50D, Hugging Face transformers,
last hidden layer). Per species and gene, the per-residue vectors are averaged over the residues that sit in the
contact columns of `mitoport burial`'s sites.tsv, and over those in the burial-matched control columns; a gap in a
column contributes nothing. Only these two means are kept (1,280 numbers each), not the per-residue vectors.

Needs a CUDA GPU (or patience on a CPU) and `uv sync --group fold`, which brings torch and transformers. The
model's weights are read from the Hugging Face cache; the commit of the checkpoint actually loaded is printed and
stored in the table's metadata. Letters ESM-2 does not know are read as X, and counted.

Writes one parquet table: slug, gene, set (contact, control), residues (how many were averaged), e0 .. e1279.
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ..longevity.aa_pgls import select
from ..longevity.structure_pgls import CLASSES
from ..seqio import read_fasta

MODEL = "facebook/esm2_t33_650M_UR50D"
COX = ("COX1", "COX2", "COX3")
ESM_LETTERS = set("LAGVSERTIDPKQNFYMHWCXBUZO")


def residue_index(row: str, cols) -> list[int]:
    """0-based positions in the ungapped sequence of `row` of the residues at alignment columns `cols` (0-based),
    gaps left out."""
    pos = np.cumsum([c != "-" for c in row]) - 1
    return [int(pos[c]) for c in cols if row[c] != "-"]


def clean(seq: str) -> tuple[str, int]:
    """The sequence with every letter ESM-2 does not know replaced by X, and how many were replaced."""
    out = "".join(c if c in ESM_LETTERS else "X" for c in seq.upper())
    return out, sum(a != b for a, b in zip(seq.upper(), out))


def species_of(data: Path, composition: Path, fourfold: Path, tree_file: Path, classes: tuple[str, ...]) -> list[str]:
    return sorted({s for cls in classes for s in select(data, composition, fourfold, tree_file, cls)[0].slug})


def main(data: Path, composition: Path, fourfold: Path, tree_file: Path, work: Path, out: Path,
         classes: tuple[str, ...] = CLASSES, model_name: str = MODEL, batch: int = 16) -> Path:
    import torch
    from transformers import AutoModel, AutoTokenizer

    sites = pd.read_csv(work / "sites.tsv", sep="\t")
    cols = {grp: {g: sites[(sites.gene == g) & (sites.group == grp)].aln_col.to_numpy() - 1 for g in COX}
            for grp in ("contact", "control")}
    slugs = species_of(data, composition, fourfold, tree_file, classes)
    alns = {g: read_fasta(data / "alignments" / f"{g}.aa.fasta.gz") for g in COX}
    device = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name, add_pooling_layer=False).to(device).eval()
    commit = getattr(model.config, "_commit_hash", None)
    print(f"{model_name} at {commit} on {device}; torch {torch.__version__}; {len(slugs)} species x {len(COX)} genes")
    jobs, replaced = [], 0
    for g in COX:
        for s in slugs:
            if s not in alns[g]:
                continue
            row = alns[g][s]; seq, k = clean(row.replace("-", "")); replaced += k
            jobs.append((s, g, seq, {grp: residue_index(row, cols[grp][g]) for grp in cols}))
    print(f"{len(jobs)} sequences; letters read as X: {replaced}")
    jobs.sort(key=lambda j: len(j[2]))
    rows, t0 = [], time.time()
    with torch.no_grad():
        for i in range(0, len(jobs), batch):
            part = jobs[i:i + batch]
            enc = tok([j[2] for j in part], return_tensors="pt", padding=True).to(device)
            h = model(**enc).last_hidden_state.float().cpu().numpy()
            for k, (s, g, seq, idx) in enumerate(part):
                per = h[k, 1:1 + len(seq)]               # drop <cls>; <eos> and padding lie beyond the sequence
                for grp, ix in idx.items():
                    v = per[ix].mean(axis=0) if ix else np.full(per.shape[1], np.nan, dtype=np.float32)
                    rows.append({"slug": s, "gene": g, "set": grp, "residues": len(ix), **{f"e{d}": v[d] for d in range(len(v))}})
    print(f"embedded in {time.time() - t0:.0f} s")
    res = pd.DataFrame(rows).sort_values(["gene", "set", "slug"]).reset_index(drop=True)
    float_cols = [c for c in res.columns if c.startswith("e") and c[1:].isdigit()]
    res[float_cols] = res[float_cols].astype(np.float32)
    out.parent.mkdir(parents=True, exist_ok=True)
    import pyarrow as pa
    import pyarrow.parquet as pq
    table = pa.Table.from_pandas(res, preserve_index=False)
    meta = {"model": model_name, "commit": str(commit), "torch": torch.__version__, "device": device}
    table = table.replace_schema_metadata({**(table.schema.metadata or {}), b"mitoage": json.dumps(meta).encode()})
    pq.write_table(table, out)
    print(f"wrote {out}: {len(res)} rows")
    return out
