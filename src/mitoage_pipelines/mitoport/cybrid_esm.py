"""H17, the language-model scores: does ESM-2 find the donor's residues unlikely in the host's protein, and more so
with the host's nuclear partner beside it?

Two readings of ESM-2 650M (facebook/esm2_t33_650M_UR50D, masked marginals), both in the host's own protein:

- alone. The host's mtDNA-encoded protein with one position masked; LLR = log p(donor's residue) - log p(host's
  residue) there. `esm_contact`, `esm_control` and `esm_all` are minus the sum of LLR over the positions where donor
  and host differ, at the contact residues, at their burial-matched controls, and over the whole proteins of the
  complex: high when the model finds the donor's residues unlikely. This is the usual zero-shot variant score; it
  never sees a nuclear subunit, so it measures how unusual the donor's protein is, not a mismatch.
- with the partner. The host's nuclear subunit, a linker of 25 glycines, then the host's mitochondrial protein, the
  same position masked. For a contact residue the LLR is averaged over the nuclear subunits it touches.
  `esm_partner_effect` is (minus the sum with the partner) minus (minus the sum alone) over the contact residues:
  what the partner's presence adds. `esm_partner_effect_control` is the same at each contact's control residue with
  the same partners, and `esm_partner_excess`, the score, is their difference. A pair longer than the model's 1,022
  residues is not assessed and is counted (`partner_pairs_too_long`).

The contact residues, their controls and their partners are those of the human structures (`mitoport burial`,
`partners.pairs`), carried to the host by the alignment columns; for a mouse host that is the human map on mouse
sequences. Host and donor proteins are H16's (`download cybrid`); the host's nuclear proteins are its rows of the
vertebrate partner alignments, gaps removed.

Needs a CUDA GPU and `uv sync --group fold`. Writes esm_llr.tsv (one row per host, gene, position, context and
donor residue) and esm_scores.tsv (one row per host, donor and complex). The scoring is separate from the model
(`score`), so that it is tested without one.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from ..seqio import read_fasta
from . import complexes, cybrid
from .burial import column_of_position
from .cybrid_coev import AA, IDX, complex_pairs

MODEL = "facebook/esm2_t33_650M_UR50D"
LINKER = "G" * 25
MAX_RESIDUES = 1022
ALONE = "alone"


def index_on_columns(human_row: str, protein: str) -> dict[int, int]:
    """{0-based alignment column: 0-based index in `protein`} for the residues of `protein` that a global alignment
    places on a residue of the human row."""
    col = column_of_position(human_row)
    human = human_row.replace("-", "").upper()
    a = cybrid.aligner().align(human, protein)[0]
    return {col[h0 + k + 1]: o0 + k for (h0, h1), (o0, o1) in zip(*a.aligned) for k in range(h1 - h0)}


def tasks(host_seq: dict[str, str], host_idx: dict[str, dict[int, int]], donors_mt: dict[str, dict[str, dict[int, str]]],
          pr: pd.DataFrame, host_nuc: dict[str, str]) -> tuple[dict, int]:
    """({(gene, context): (sequence, offset of the mitochondrial protein in it, sorted protein indices to mask)},
    pairs too long). Context is ALONE or a nuclear gene. A position is masked when at least one donor carries an
    amino acid there that differs from the host's."""
    differs = {}
    for g, seq in host_seq.items():
        cols = set()
        for d in donors_mt.values():
            for c, aa in d[g].items():
                i = host_idx[g].get(c)
                if i is not None and aa in IDX and seq[i] in IDX and aa != seq[i]:
                    cols.add(c)
        differs[g] = cols
    out = {(g, ALONE): (seq, 0, sorted(host_idx[g][c] for c in differs[g])) for g, seq in host_seq.items() if differs[g]}
    too_long = 0
    for (g, n), grp in pr.groupby(["mt_gene", "nuc_gene"], sort=True):
        cols = (set(grp.mt_col) | set(grp.ctl_col)) & differs[g]
        if not cols or n not in host_nuc:
            continue
        seq = host_nuc[n] + LINKER + host_seq[g]
        if len(seq) > MAX_RESIDUES:
            too_long += 1
            continue
        out[(g, n)] = (seq, len(host_nuc[n]) + len(LINKER), sorted(host_idx[g][c] for c in cols))
    return out, too_long


def score(pr: pd.DataFrame, contact_cols: dict[str, set[int]], control_cols: dict[str, set[int]],
          host_seq: dict[str, str], host_idx: dict[str, dict[int, int]], donor_mt: dict[str, dict[int, str]],
          logp: dict[tuple[str, str, int], np.ndarray]) -> dict:
    """The sums for one (host, donor) on one complex. `logp[(gene, context, protein index)]` is the model's 20
    log-probabilities (in the order of AA) at that masked position."""
    def llr(g, ctx, c):
        """LLR of the donor's residue at column c; 0 where it is the host's; None where it cannot be read."""
        i = host_idx[g].get(c)
        if i is None:
            return None
        h, d = host_seq[g][i], donor_mt[g].get(c, "?")
        if d not in IDX or h not in IDX:
            return None
        if d == h:
            return 0.0
        v = logp.get((g, ctx, i))
        return None if v is None else float(v[IDX[d]] - v[IDX[h]])

    out = dict(esm_contact=0.0, esm_control=0.0, esm_all=0.0, positions_unknown=0, n_contact=0)
    for g in host_seq:
        for c in sorted(host_idx[g]):
            v = llr(g, ALONE, c)
            if v is None:
                out["positions_unknown"] += 1
                continue
            out["esm_all"] -= v
            if c in contact_cols[g]:
                out["esm_contact"] -= v; out["n_contact"] += 1
            if c in control_cols[g]:
                out["esm_control"] -= v
    out["esm_contact_mean"] = out["esm_contact"] / out["n_contact"] if out["n_contact"] else np.nan
    eff = {"mt_col": 0.0, "ctl_col": 0.0}; assessed = 0
    for (g, c, ctl), grp in pr.groupby(["mt_gene", "mt_col", "ctl_col"], sort=True):
        vals = {}
        for name, col in (("mt_col", c), ("ctl_col", ctl)):
            alone = llr(g, ALONE, col)
            with_p = [llr(g, n, col) for n in sorted(set(grp.nuc_gene))]
            with_p = [v for v in with_p if v is not None]
            if alone is None or not with_p:
                vals = None; break
            vals[name] = -(float(np.mean(with_p)) - alone)
        if vals is None:
            continue
        assessed += 1
        for k in eff:
            eff[k] += vals[k]
    out["partner_residues_assessed"] = assessed
    out["esm_partner_effect"] = eff["mt_col"] if assessed else np.nan
    out["esm_partner_effect_control"] = eff["ctl_col"] if assessed else np.nan
    out["esm_partner_excess"] = eff["mt_col"] - eff["ctl_col"] if assessed else np.nan
    return out


def model_logprobs(model_name: str = MODEL, batch: int = 8):
    """A function (sequence, indices) -> array [len(indices), 20] of log-probabilities at each index masked in turn,
    and the commit of the checkpoint loaded."""
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    tok = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForMaskedLM.from_pretrained(model_name).to(device).eval()
    aa_ids = torch.tensor(tok.convert_tokens_to_ids(list(AA)), device=device)
    commit = getattr(model.config, "_commit_hash", None) or "unknown"
    print(f"{model_name} at {commit} on {device}")

    def run(seq: str, idx: list[int]) -> np.ndarray:
        ids = tok(seq, return_tensors="pt")["input_ids"].to(device)
        assert ids.shape[1] == len(seq) + 2, "one token per residue plus the two ends"
        out = []
        with torch.no_grad():
            for k in range(0, len(idx), batch):
                part = idx[k:k + batch]
                x = ids.repeat(len(part), 1)
                rows = torch.arange(len(part), device=device); cols = torch.tensor(part, device=device) + 1
                x[rows, cols] = tok.mask_token_id
                logits = model(input_ids=x).logits[rows, cols]
                out.append(torch.log_softmax(logits.float(), dim=-1)[:, aa_ids].cpu().numpy())
        return np.concatenate(out) if out else np.zeros((0, 20))
    return run, commit


def main(data: Path, gb: Path, nucaln: Path, out: Path, measurements: Path = cybrid.MEASUREMENTS,
         logprobs=None) -> pd.DataFrame:
    out.mkdir(parents=True, exist_ok=True)
    meas = cybrid.read_measurements(measurements)
    records = cybrid.read_proteins(gb)
    species = sorted(set(meas.host) | set(meas.donor))
    prot = {s: cybrid.species_proteins(records, s) for s in species}
    commit = "injected"
    if logprobs is None:
        logprobs, commit = model_logprobs()
    from .cybrid_coev import on_columns
    rows, llr_rows = [], []
    for name, cx in complexes.COMPLEXES.items():
        pr, nuc_alns = complex_pairs(cx, nucaln)
        sites = pd.read_csv(cx.work / "sites.tsv", sep="\t", dtype={"chain": str})
        cols = {grp: {g: set(sites.aln_col[(sites.gene == g) & (sites.group == grp)] - 1) for g in cx.mt}
                for grp in ("contact", "control")}
        human_row = {g: read_fasta(data / "alignments" / f"{g}.aa.fasta.gz")["homo-sapiens"] for g in cx.mt}
        for h in sorted(set(meas.host)):
            if any(prot[h][g][1] for g in cx.mt):
                raise SystemExit(f"{h}: a host protein of complex {name} is a partial record")
            host_seq = {g: prot[h][g][0] for g in cx.mt}
            host_idx = {g: index_on_columns(human_row[g], host_seq[g]) for g in cx.mt}
            host_nuc = {n: a[h].replace("-", "").upper() for n, a in nuc_alns.items() if h in a}
            donors = sorted(set(meas.donor[meas.host == h]))
            donors_mt = {d: {g: on_columns(human_row[g], prot[d][g][0], prot[d][g][1]) for g in cx.mt} for d in donors}
            todo, too_long = tasks(host_seq, host_idx, donors_mt, pr, host_nuc)
            print(f"complex {name}, host {h}: {len(todo)} sequences, {sum(len(t[2]) for t in todo.values())} masked "
                  f"positions, {too_long} partner pairs over {MAX_RESIDUES} residues (not assessed)")
            logp = {}
            for (g, ctx), (seq, off, idx) in sorted(todo.items()):
                lp = logprobs(seq, [off + i for i in idx])
                for i, v in zip(idx, lp):
                    logp[(g, ctx, i)] = v
            for d in donors:
                rows.append({"host": h, "donor": d, "complex": name, "partner_pairs_too_long": too_long,
                             **score(pr, cols["contact"], cols["control"], host_seq, host_idx, donors_mt[d], logp)})
                for g in cx.mt:
                    for c, i in sorted(host_idx[g].items()):
                        aa = donors_mt[d][g].get(c, "?")
                        if aa in IDX and host_seq[g][i] in IDX and aa != host_seq[g][i]:
                            for ctx in sorted({k[1] for k in logp if k[0] == g and k[2] == i}):
                                v = logp[(g, ctx, i)]
                                llr_rows.append({"host": h, "donor": d, "complex": name, "gene": g, "aln_col": c + 1,
                                                 "host_index": i + 1, "host_aa": host_seq[g][i], "donor_aa": aa, "context": ctx,
                                                 "llr": float(v[IDX[aa]] - v[IDX[host_seq[g][i]]])})
    res = pd.DataFrame(rows)
    res.attrs["checkpoint"] = commit
    res.to_csv(out / "esm_scores.tsv", sep="\t", index=False, float_format="%.5f")
    pd.DataFrame(llr_rows).to_csv(out / "esm_llr.tsv", sep="\t", index=False, float_format="%.5f")
    (out / "esm_checkpoint.txt").write_text(f"{MODEL}\t{commit}\n")
    with pd.option_context("display.width", 220, "display.max_rows", 200):
        print(res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("wrote", out / "esm_scores.tsv", out / "esm_llr.tsv")
    return res
