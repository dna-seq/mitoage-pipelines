"""Predicted structures of the 13 mitochondrially encoded proteins for a set of species, with ESMFold2 (Biohub;
the single-sequence checkpoint biohub/ESMFold2-Fast, one diffusion sample per protein).

Needs a CUDA GPU and the `fold` dependency group (`uv sync --group fold`: esm, torch). One structure per gene
and species at <out>/<GENE>/<species>.pdb; an existing file is kept, so a run resumes after a time limit.
Sequences are the ungapped rows of the MitoAge alignments, X removed. The mean pLDDT of each model is appended
to <out>/plddt.tsv (no header: gene, species, length, mean pLDDT).

The published run (2026-11-15, one A100 40 GB, esm 3.4.3, torch 2.15.0, Python 3.14) folded the 62 mammals of
`mitoport align`'s species.txt in 11,949 s. The esm package exposes several call paths and only some write a
PDB from half-precision tensors, so the one that works is settled on a short sequence first and printed.
Ported from the cluster script `struct/fold.py` without change of method.
"""
import time
import traceback
from pathlib import Path

import numpy as np

from ..genes import MT_GENES
from ..seqio import read_fasta

CHECKPOINT = "biohub/ESMFold2-Fast"
GENES = "ND1 ND2 ND3 ND4 ND4L ND5 ND6 CYTB COX1 COX2 COX3 ATP6 ATP8".split()     # the published order
assert sorted(GENES) == sorted(MT_GENES)


def mean_plddt(pdb_text: str) -> float:
    """Mean of the B-factor column (pLDDT) over the CA atoms of a PDB text."""
    v = [float(line[60:66]) for line in pdb_text.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CA"]
    return float(np.mean(v)) if v else float("nan")


def main(build: Path, species_file: Path, out: Path, checkpoint: str = CHECKPOINT, genes: list[str] | None = None,
         samples: int = 1) -> None:
    import torch
    import esm.models.esmfold2 as F

    genes = genes or GENES
    out.mkdir(parents=True, exist_ok=True)
    kw = dict(num_diffusion_samples=samples)     # the checkpoint's default is 8: eight times memory and time
    sp = [s for s in species_file.read_text().split() if s]
    print(f"{len(sp)} species x {len(genes)} genes, checkpoint {checkpoint}", flush=True)
    if not torch.cuda.is_available():
        raise RuntimeError("ESMFold2 needs a CUDA GPU; none is visible")
    free, total = torch.cuda.mem_get_info()
    dt = torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16
    print(f"card: {torch.cuda.get_device_name(0)}, {free / 2**30:.1f} GiB free of {total / 2**30:.1f}; weights in {dt}", flush=True)
    model = None
    for ctor in (F.EsmFold2Model, F.EsmFold2HFAdapter):
        try:
            model = ctor.from_pretrained(checkpoint).cuda().eval(); print(f"loaded through {ctor.__name__}", flush=True); break
        except Exception:
            print(f"{ctor.__name__}.from_pretrained failed:"); traceback.print_exc()
    if model is None:
        raise SystemExit("no ESMFold2 checkpoint could be loaded")
    try:
        model = model.to(dt)
    except Exception as e:
        print(f"half precision refused ({type(e).__name__}), staying in fp32", flush=True); dt = torch.float32
    try:
        model.set_chunk_size(64)
    except Exception as e:
        print(f"chunking not available ({type(e).__name__})", flush=True)

    def path_pdb(seq):
        return model.infer_protein_as_pdb(seq, **kw)

    def path_input(seq):          # the writer cannot read half-precision tensors
        o = model.infer_protein(seq, **kw)
        o = {k: (v.float() if torch.is_tensor(v) and v.dtype in (torch.bfloat16, torch.float16) else v) for k, v in o.items()}
        return model.output_to_pdb(o)
    call = None
    for name, fn in [("infer_protein + fp32 output", path_input), ("infer_protein_as_pdb", path_pdb)]:
        try:
            with torch.no_grad(), torch.autocast("cuda", dtype=dt, enabled=dt != torch.float32):
                s = fn("MKTFFVLLLAIGSLAQSVLAGKVFGRCELA")
            if isinstance(s, str) and len(s) > 200:
                call = fn; print(f"call path: {name}", flush=True); break
            print(f"call path {name} returned {type(s)}, not a structure text")
        except Exception:
            print(f"call path {name} failed:"); traceback.print_exc()
    if call is None:
        raise SystemExit("no working ESMFold2 call path")

    done = fail = skip = 0; t0 = time.time()
    with open(out / "plddt.tsv", "a") as log:
        for g in genes:
            seqs = {k: v.replace("-", "").replace("X", "") for k, v in read_fasta(build / "alignments" / f"{g}.aa.fasta.gz").items()}
            (out / g).mkdir(exist_ok=True); tg = time.time()
            for s in sp:
                p = out / g / f"{s}.pdb"
                if p.exists() and p.stat().st_size > 1000:
                    skip += 1; continue
                if s not in seqs or len(seqs[s]) < 30:
                    continue
                try:
                    with torch.no_grad(), torch.autocast("cuda", dtype=dt, enabled=dt != torch.float32):
                        text = call(seqs[s])
                    p.write_text(text); done += 1
                    log.write(f"{g}\t{s}\t{len(seqs[s])}\t{mean_plddt(text):.2f}\n"); log.flush()
                except Exception as e:
                    print(f"  failed {g} {s} ({len(seqs[s])} aa): {type(e).__name__} {e}", flush=True); fail += 1
                torch.cuda.empty_cache()
            print(f"{g}: {len(list((out / g).glob('*.pdb')))} structures ({time.time() - tg:.0f} s for this gene, "
                  f"{time.time() - t0:.0f} s total)", flush=True)
    print(f"folded {done}, kept {skip}, failed {fail}, {time.time() - t0:.0f} s")
    if fail:
        raise SystemExit(f"{fail} structures failed; rerun to retry them (finished ones are kept)")
