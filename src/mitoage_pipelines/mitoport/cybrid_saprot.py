"""H36, step 3: a protein language model that sees a mitochondrial chain and one nuclear chain together. Which nuclear subunit's
ortholog restores the fit of the donor's mitochondrial chain, and which nuclear residues carry it?

The model is SaProt 650M (`westlake-repl/SaProt_650M_PDB`, MIT): ESM-2's architecture trained on residues paired with the
3Di structure token of Foldseek, here the weights tuned on solved structures. Its input is one token per residue, the amino
acid and its structure token ("Mp": methionine, structure state p); a masked residue is "#p", its structure token kept. The
structure tokens come from our own solved complexes by Foldseek's binary at a pinned release (`structureto3didescriptor`).
Nothing here sees an alignment, a species or a genetic code; it sees the host's backbone as solved.

What is read, fixed in docs/HYPOTHESES.md (H36, "Step 3 made exact", 2026-11-25) before any score:

- One input is one mitochondrial chain, a linker of LINKER glycines with no structure ("G#"), and one nuclear chain, in the
  system of the swap (`cybrid_swap`: the complex's mtDNA-encoded chains and the chains of the nuclear gene u that touch them,
  on the host's own solved structure). The pairs are the (mitochondrial chain, chain of u) that have a heavy atom of each
  within CONTACT_A. SaProt was trained on single chains, so a joined pair is outside what it was built for: what it does with
  the second chain is tested by shuffling it (`shuffled`), not assumed. A pair longer than MAX_RESIDUES residues is not
  assessed and counted; a chain whose structure tokens do not match its sequence is not assessed.
- The score of an input is the pseudo-log-likelihood, in nats, of the mitochondrial chain's residues that touch u: each such
  residue in turn is masked ("#" and its own structure token), and its log-probability is that of the residue the input holds,
  taken over the 20 amino acids at the structure token it has (a softmax over those 20 logits). The sum over the residues and
  over the pairs of u is the task's PLL. The model is deterministic: no seed, no sweep.
- The mitochondrial chain carries the donor's substitutions (the donor read through the host, as H28 builds them). B0: beside the
  host's own u; B1: beside u carrying the donor's ortholog; restoration = PLL(B1) - PLL(B0) (positive: the donor's own u makes the
  donor's mitochondrial residues likelier). The tables keep the step-2 rule's own code by storing the misfit, minus the PLL, in the
  column the rule reads, so that "B0 minus B1" is the restoration there as here.
- The null: the same B1 with the ortholog of a random other mammal of the alignment, as the swap draws it (within 25% of the
  donor's differences at the contact residues, at least 3 else not assessed, 50 draws, seeds (0, draw)); not run when the donor's
  ortholog changes no placed residue (the restoration is exactly 0).
- Which residues matter: for every residue of u that the ortholog changes, that residue alone changed (host to donor) beside the
  donor's mitochondrial chain, the change of the PLL (`residues`); beside it the residue's distance to the mitochondrial chains
  and whether it is at the interface, so that a model that does not see the interface shows it. Attention is not read.
- The partner is shuffled SHUFFLES times (the residues of u permuted among its amino-acid places, structure tokens kept, seeds
  (1, draw)): PLL(B0) minus the shuffled PLL says how much the model uses the real partner.

`saprot` runs one (donor, complex, gene) and writes saprot/<donor>.<complex>.<gene>.tsv (one row per kind and draw),
.pairs.tsv (one row per pair and read), .residues.tsv and .cases.tsv, the last one written last: a task whose case table is there
is not run again (a requeued job resumes). `saprot-collect` joins them into the tables and applies the step-2 rule.
"""
import hashlib
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import logsumexp

from . import cybrid, cybrid_stab as cs, cybrid_swap as sw
from .cybrid_energy import STRUCTURES, THREE, insertions
from .complexes import COMPLEXES

MODEL_REPO = "westlake-repl/SaProt_650M_PDB"
MODEL_REVISION = "6f7468fdd0cdffe7e780bd4f57fa291df6b01867"
MODEL_FILES = {         # file -> SHA-256
    "pytorch_model.bin": "f37c51299a3ee31ead2823b835300e5d12f60453a2f69d2b977aa2a36b99641a",
    "config.json": "927bf8961c93e0667cecb6f928aa7077e1a7f55c8483ae58173f155261f1c3b3",
    "vocab.txt": "d4289b3ff55d50b1aa33a3141ec344fc2a69cb3dca43aeefe9fccfee371ffd28",
}
FOLDSEEK_RELEASE = "10-941cd33"
FOLDSEEK_URL = f"https://github.com/steineggerlab/foldseek/releases/download/{FOLDSEEK_RELEASE}/foldseek-linux-avx2.tar.gz"
FOLDSEEK_SHA256 = "e6a3e3860e76b0b7b858736f627fe680d4d7f0fc7298301d8a92b9f711d74ac7"
FOLDSEEK_VERSION = "941cd33ff0771cd2e3f144e3293e22a2b87e9fda"

AA = "ACDEFGHIKLMNPQRSTVWY"
LINKER = 25
MAX_RESIDUES = 1022                 # the model's training length of 1,024 tokens, less its two end tokens
BATCH = 16                          # masked copies of one input read in one pass
SHUFFLES = 3
SEED = 0                            # the model is deterministic: one "seed" in the tables, so that the swap's tables read them
SCORE_COLUMNS = ["host", "donor", "complex", "gene", "kind", "draw", "species", "seed", "substituted", "pairs", "positions", "pll", "ddG"]
PAIR_COLUMNS = ["host", "donor", "complex", "gene", "kind", "draw", "mt_gene", "mt_chain", "u_chain", "positions", "pll"]
RESIDUE_COLUMNS = ["host", "donor", "complex", "gene", "mt_gene", "mt_chain", "u_chain", "resnum", "host_aa", "donor_aa", "interface",
                   "min_dist", "positions", "delta_pll", "top_mt_resnum", "top_mt_delta"]
NUMERIC_SCORES = ["draw", "seed", "substituted", "pairs", "positions", "pll", "ddG"]
NUMERIC_PAIRS = ["draw", "positions", "pll"]
NUMERIC_RESIDUES = ["resnum", "min_dist", "positions", "delta_pll", "top_mt_resnum", "top_mt_delta"]
STATE = "saprot"                    # the tasks' folder under the output folder


# --- the pinned downloads -------------------------------------------------------------------------------------------


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def download_weights(dest: Path, manifest: Path) -> list[Path]:
    """SaProt 650M (PDB) at the pinned revision into `dest`; a file whose SHA-256 is not the pinned one is refused and removed."""
    from huggingface_hub import hf_hub_download
    from ..download import append_manifest
    dest.mkdir(parents=True, exist_ok=True); got = []
    for name, digest in MODEL_FILES.items():
        p = dest / name
        if not (p.is_file() and sha256(p) == digest):
            p = Path(hf_hub_download(MODEL_REPO, name, revision=MODEL_REVISION, local_dir=dest))
            if sha256(p) != digest:
                p.unlink()
                raise SystemExit(f"{name} at {MODEL_REVISION[:7]} is not the pinned file (SHA-256 differs): removed")
        got.append(p); print(f"{name}  sha256 {digest}")
    append_manifest(manifest, "huggingface", f"{MODEL_REPO}@{MODEL_REVISION}", dest, sha256=MODEL_FILES)
    return got


def download_foldseek(dest: Path, manifest: Path) -> Path:
    """Foldseek's static binary at the pinned release into `dest`/foldseek; the archive is kept only with its pinned SHA-256 and the
    binary only if it reports the pinned build."""
    import tarfile
    from ..download import _request, append_manifest
    dest.mkdir(parents=True, exist_ok=True)
    binary = dest / "foldseek"
    if not (binary.is_file() and foldseek_version(binary) == FOLDSEEK_VERSION):
        body = _request(FOLDSEEK_URL, accept="*/*")
        if hashlib.sha256(body).hexdigest() != FOLDSEEK_SHA256:
            raise SystemExit(f"the archive at {FOLDSEEK_URL} is not the pinned file (SHA-256 differs): nothing written")
        tar = dest / "foldseek-linux-avx2.tar.gz"
        tar.write_bytes(body)
        with tarfile.open(tar) as t:
            member = t.getmember("foldseek/bin/foldseek")
            binary.write_bytes(t.extractfile(member).read())
        binary.chmod(0o755)
        if foldseek_version(binary) != FOLDSEEK_VERSION:
            binary.unlink()
            raise SystemExit("the binary does not report the pinned build")
    append_manifest(manifest, "github", f"steineggerlab/foldseek@{FOLDSEEK_RELEASE}", dest, sha256={"foldseek-linux-avx2.tar.gz": FOLDSEEK_SHA256})
    return binary


def foldseek_version(binary: Path) -> str:
    return subprocess.run([str(binary), "version"], capture_output=True, text=True, check=True).stdout.strip()


# --- structure tokens -----------------------------------------------------------------------------------------------


def parse_foldseek(text: str) -> dict[str, tuple[str, str]]:
    """{chain id: (amino acids, 3Di states in lower case)} from `structureto3didescriptor --chain-name-mode 1` output: one line per
    chain, the first field "<file>_<chain> <title>", the second the residues, the third the states. A chain seen twice is refused."""
    out: dict[str, tuple[str, str]] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 3:
            raise SystemExit(f"a line of Foldseek's output has {len(f)} fields, not 3: {line[:60]!r}")
        chain = f[0].split(" ")[0].rsplit("_", 1)[1]
        if chain in out:
            raise SystemExit(f"chain {chain} appears twice in Foldseek's output")
        if len(f[1]) != len(f[2]):
            raise SystemExit(f"chain {chain}: {len(f[1])} residues and {len(f[2])} structure states")
        out[chain] = (f[1], f[2].lower())
    return out


def structure_tokens(cif: Path, binary: Path | None, cache: Path) -> dict[str, tuple[str, str]]:
    """The tokens of every chain of `cif`, from the cache table when it is there, else from Foldseek's binary (written to the cache
    as `<stem>.3di.tsv`: chain, residues, states)."""
    f = cache / f"{cif.stem}.3di.tsv"
    if f.is_file():
        t = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        return {r.chain: (r.aa, r.di) for r in t.itertuples()}
    if binary is None or not binary.is_file():
        raise SystemExit(f"no tokens of {cif.stem} in {cache} and no Foldseek binary: `mitoage download foldseek`, `mitoage mitoport saprot-tokens`")
    cache.mkdir(parents=True, exist_ok=True)
    work = cache / f".{cif.stem}.work"
    work.mkdir(exist_ok=True)
    subprocess.run([str(binary), "structureto3didescriptor", "-v", "0", "--threads", "1", "--chain-name-mode", "1", str(cif), str(work / "out")],
                   check=True, capture_output=True, text=True, cwd=work)
    tok = parse_foldseek((work / "out").read_text())
    for p in work.iterdir():
        p.unlink()
    work.rmdir()
    pd.DataFrame([{"chain": k, "aa": a, "di": d} for k, (a, d) in tok.items()]).to_csv(f, sep="\t", index=False)
    return tok


def chain_states(chain_seq: str, tok_aa: str, tok_di: str) -> str | None:
    """The 3Di states of a chain as `read_chains` holds it, or None when Foldseek's chain is not that chain: the same length and the
    same residue at every place, except that Foldseek writes X where it does not know a residue (a modified one) and the chain holds
    its standard parent."""
    if len(tok_aa) != len(chain_seq) or len(tok_di) != len(chain_seq):
        return None
    if any(a != b and b != "X" for a, b in zip(chain_seq, tok_aa)):
        return None
    return tok_di


# --- the input ------------------------------------------------------------------------------------------------------


def read_vocab(path: Path) -> dict[str, int]:
    """{token: id} in the order of the file (the model's embedding rows); the specials come first."""
    v = [t for t in path.read_text().split("\n") if t]
    if len(set(v)) != len(v) or v[:5] != ["<cls>", "<pad>", "<eos>", "<unk>", "<mask>"]:
        raise SystemExit(f"{path} is not SaProt's vocabulary")
    return {t: i for i, t in enumerate(v)}


def sa_tokens(seq: str, states: str) -> list[str]:
    """One token per residue: its amino acid and its structure state; a residue that is not one of the 20 is "#" (no residue) and its state."""
    if len(seq) != len(states):
        raise ValueError(f"{len(seq)} residues and {len(states)} states")
    return [(a if a in AA else "#") + s for a, s in zip(seq, states)]


def apply(seq: str, subs: dict[int, str]) -> str:
    """`seq` with the residues at `subs` replaced; a place that holds no amino acid, or a replacement that is none, is refused."""
    s = list(seq)
    for i, aa in subs.items():
        if s[i] not in AA or aa not in AA:
            raise ValueError(f"substitution at a place that holds no amino acid: {s[i]} to {aa} at {i}")
        s[i] = aa
    return "".join(s)


def joined(m_seq: str, m_states: str, u_seq: str, u_states: str, m_subs: dict[int, str] | None = None,
           u_subs: dict[int, str] | None = None) -> tuple[list[str], int]:
    """(the tokens of the mitochondrial chain, the linker and the nuclear chain, the index of the nuclear chain's first token), the
    substitutions applied to each chain."""
    m = sa_tokens(apply(m_seq, m_subs or {}), m_states)
    u = sa_tokens(apply(u_seq, u_subs or {}), u_states)
    return m + ["G#"] * LINKER + u, len(m) + LINKER


def shuffled(u_seq: str, draw: int) -> dict[int, str]:
    """{index: residue} of the places that change when the amino acids of `u_seq` are permuted among the places that hold one,
    from the generator seeded (1, draw)."""
    idx = [i for i, a in enumerate(u_seq) if a in AA]
    perm = np.random.default_rng((1, draw)).permutation(len(idx))
    return {idx[k]: u_seq[idx[p]] for k, p in enumerate(perm) if u_seq[idx[p]] != u_seq[idx[k]]}


def masked_logp(logits_at, vocab: dict[str, int], tokens: list[str], positions: list[int], batch: int = BATCH) -> np.ndarray:
    """The log-probability of the residue each of `positions` holds, one place masked at a time, for one input.
    `logits_at(ids [B, T], cols [B])` returns the model's logits at column cols[b] of row b ([B, vocab]); the input is the tokens between
    <cls> and <eos>, a column is its position plus one. A masked place is "#" and its own structure state; the log-probability is a
    softmax over the 20 amino acids at that state. A place that holds no amino acid (a token starting "#") has no value (NaN)."""
    ids = np.array([vocab["<cls>"], *(vocab[t] for t in tokens), vocab["<eos>"]], dtype=np.int64)
    out = np.full(len(positions), np.nan)
    todo = [k for k, p in enumerate(positions) if tokens[p][0] in AA]
    for b in range(0, len(todo), batch):
        part = todo[b:b + batch]
        x = np.tile(ids, (len(part), 1))
        cols = np.array([positions[k] + 1 for k in part])
        state = [tokens[positions[k]][1] for k in part]
        x[np.arange(len(part)), cols] = [vocab["#" + s] for s in state]
        logits = np.asarray(logits_at(x, cols), dtype=np.float64)
        for j, k in enumerate(part):
            lp = logits[j, [vocab[a + state[j]] for a in AA]]
            out[k] = (lp - logsumexp(lp))[AA.index(tokens[positions[k]][0])]
    return out


def load_model(weights_dir: Path, device: str):
    """(logits_at, the revision): SaProt from `weights_dir`, each file checked against its pinned SHA-256, in 32-bit floats (the card
    has no bfloat16), on `device`. The GPU is touched at once."""
    import torch
    from transformers import EsmForMaskedLM
    for name, digest in MODEL_FILES.items():
        if not (weights_dir / name).is_file() or sha256(weights_dir / name) != digest:
            raise SystemExit(f"{weights_dir / name} is missing or is not the pinned file: `mitoage download saprot`")
    if device.startswith("cuda"):
        torch.cuda.init()
    model = EsmForMaskedLM.from_pretrained(str(weights_dir), dtype=torch.float32).to(device).eval()
    print(f"{MODEL_REPO}@{MODEL_REVISION[:7]} on {device}, attention {model.config._attn_implementation}", flush=True)

    def logits_at(x: np.ndarray, cols: np.ndarray) -> np.ndarray:
        with torch.no_grad():
            t = torch.from_numpy(x).to(device)
            out = model(input_ids=t, attention_mask=torch.ones_like(t)).logits
            return out[torch.arange(len(cols), device=device), torch.from_numpy(cols).to(device)].float().cpu().numpy()
    return logits_at, MODEL_REVISION


# --- one task -------------------------------------------------------------------------------------------------------


def residue_distances(chains: list[dict], mt_idx: list[int], c: int) -> np.ndarray:
    """The distance of each residue of chain `c` to the nearest heavy atom of the mtDNA-encoded chains `mt_idx`."""
    from scipy.spatial import cKDTree
    d = cKDTree(np.concatenate([chains[i]["heavy"] for i in mt_idx])).query(chains[c]["heavy"])[0]
    edges = np.concatenate([[0], np.cumsum(chains[c]["n_heavy"])])
    return np.array([d[a:b].min() for a, b in zip(edges[:-1], edges[1:])])


def pair_positions(chains: list[dict], mt_idx: list[int], c: int) -> dict[int, list[int]]:
    """{mitochondrial chain: its residues with a heavy atom within CONTACT_A of chain c}, for the chains that have one: the pairs of c."""
    out = {}
    for k in mt_idx:
        flags = sw.contact_flags(chains, [c], [k])[k]
        if flags.any():
            out[k] = [i for i, f in enumerate(flags) if f]
    return out


def chain_state_table(chains: list[dict], tokens: dict[str, tuple[str, str]]) -> dict[int, str | None]:
    """{chain index: its 3Di states, or None where Foldseek's chain is not the chain}."""
    return {i: chain_states(c["seq"], *tokens[c["id"]]) if c["id"] in tokens else None for i, c in enumerate(chains)}


def run(donor: str, complex_: str, gene: str, gb: Path, pdb_dir: Path, aln_dir: Path, species: Path, tokens_dir: Path, weights_dir: Path,
        out: Path, binary: Path | None = None, device: str = "cuda", draws: int = sw.DRAWS, batch: int = BATCH, shuffles: int = SHUFFLES,
        logits_for=None, vocab_for=None, donor_source=None, host: str = sw.HOST) -> tuple[pd.DataFrame, dict] | None:
    """One (donor, complex, nuclear gene): writes the task's tables under `out`/saprot and returns (the scores, the case row), or None when
    its case table is already there (skipped: a requeued job resumes). The preamble, up to the sets of substitutions, is the swap's
    (`cybrid_swap.run`): the same host structure, the same chains of u, the same placement of the orthologs, the same null candidates.
    `logits_for()` replaces the model, `vocab_for()` the vocabulary (for tests)."""
    if host != sw.HOST and donor_source is None:
        raise SystemExit(f"a host other than {sw.HOST} reads its subunits from the pinned RefSeq orthologs: give the donor orthologs")
    pid = STRUCTURES.get(host, {}).get(complex_)
    d = out / STATE; d.mkdir(parents=True, exist_ok=True)
    stem = d / f"{donor}.{complex_}.{gene}"
    if Path(f"{stem}.cases.tsv").is_file():
        print("skipped (done):", stem)
        return None

    def done(case: dict, scores: pd.DataFrame | None = None, pairs: pd.DataFrame | None = None, residues: pd.DataFrame | None = None):
        scores = pd.DataFrame(columns=SCORE_COLUMNS) if scores is None else scores
        for name, t, cols in (("", scores, SCORE_COLUMNS), (".pairs", pairs, PAIR_COLUMNS), (".residues", residues, RESIDUE_COLUMNS)):
            t = pd.DataFrame(columns=cols) if t is None else t
            t.to_csv(f"{stem}{name}.tmp", sep="\t", index=False, float_format="%.6f")
            Path(f"{stem}{name}.tmp").replace(f"{stem}{name}.tsv")
        pd.DataFrame([case]).to_csv(f"{stem}.cases.tmp", sep="\t", index=False, float_format="%.6f")
        Path(f"{stem}.cases.tmp").replace(f"{stem}.cases.tsv")          # the marker, last
        print("wrote", f"{stem}.tsv", case["status"], flush=True)
        return scores, case

    if pid is None:
        return done(sw.case_row(donor, complex_, gene, f"not assessed: no solved structure of the {host} complex {complex_}", host=host))
    cif = pdb_dir / f"{pid}.cif"
    if not cif.is_file():
        raise SystemExit(f"{cif} is missing: `mitoage download pdb {pid}`")
    aln_file = aln_dir / f"{gene}.afa"
    host_ortholog = None
    if host != sw.HOST:
        host_ortholog, reason = donor_source.get(host, gene)
        if host_ortholog is None:
            why = "no pinned source for the host" if reason is None else reason.replace("not assessed: ", "").replace("of the donor", "of the host")
            return done(sw.case_row(donor, complex_, gene, f"not assessed: {why}", host=host))
    ortholog, reason = donor_source.get(donor, gene) if donor_source is not None else (None, None)
    if reason is not None:
        return done(sw.case_row(donor, complex_, gene, reason, host=host))
    if ortholog is None and host != sw.HOST:
        return done(sw.case_row(donor, complex_, gene, "not assessed: the pinned table does not cover the donor", host=host))
    rows = sw.read_alignment(aln_file) if aln_file.is_file() else None
    if rows is not None and host not in rows:
        raise SystemExit(f"{aln_file} holds no row of {host}")
    if ortholog is None:
        if rows is None:
            return done(sw.case_row(donor, complex_, gene, "not assessed: no alignment of the gene"))
        if donor not in rows:
            return done(sw.case_row(donor, complex_, gene, "not assessed: no ortholog of the donor in the alignment"))
    records = cybrid.read_proteins(gb)
    hp, dp = (cybrid.species_proteins(records, s) for s in (host, donor))
    genes = tuple(COMPLEXES[complex_].mt)
    mt = {g: hp[g][0] for g in genes}
    chains, mt_chain = cs.system(cs.read_chains(cif), mt, {g: p for g, (p, _, _) in hp.items() if g not in genes})
    if host == sw.HOST:
        ct = pd.read_csv(COMPLEXES[complex_].work / "chains.tsv", sep="\t", dtype={"chain": str})
        by_gene = sw.gene_chains(chains, mt_chain, dict(zip(ct.chain, ct.gene)))
        no_chain = "not assessed: no chain of the gene touches the mitochondrial chains"
    else:
        nuclear = {g: o.protein for g in COMPLEXES[complex_].nuclear for o in [donor_source.get(host, g)[0]] if o is not None}
        by_gene = sw.gene_chains_by_sequence(chains, mt_chain, nuclear)
        no_chain = (f"not assessed: no chain of {pid} touching the mitochondrial chains is the host's {gene} (the structure holds none "
                    f"that is 90% identical to its RefSeq protein: it touches none, or it holds another isoform)")
    if gene not in by_gene:
        return done(sw.case_row(donor, complex_, gene, no_chain, host=host))
    u_idx, one_idx = by_gene[gene], sorted(mt_chain)
    extra: dict = {}
    if host != sw.HOST:
        differs = sum(len(sw.protein_subs(chains[c]["seq"], host_ortholog.protein)[0]) for c in u_idx)
        if differs:
            return done(sw.case_row(donor, complex_, gene, f"not assessed: the chains of {pid} differ from the host's RefSeq protein at {differs} placed residues, "
                                    "so the chain as solved is not the host's own subunit", host=host))
        extra = {"host_protein": host_ortholog.accession, "host_differs": differs}
    flags = sw.contact_flags(chains, one_idx, u_idx)
    cols = {c: sw.chain_columns(chains[c]["seq"], rows[host]) for c in u_idx} if rows is not None else None
    columns = sw.contact_columns(cols, flags) if cols is not None else []
    sp = pd.read_parquet(species, columns=["slug", "class"])
    pool = set(sp.slug[sp["class"] == sw.CLASS])

    def u_subs(row: str) -> tuple[dict[tuple[int, int], str], int, int]:
        sub, placed, lost = {}, 0, 0
        for c in u_idx:
            s, p, l = sw.ortholog_subs(chains[c]["seq"], cols[c], row)
            sub.update({(c, i): aa for i, aa in s.items()}); placed += p; lost += l
        return sub, placed, lost

    if ortholog is None:
        donor_sub, placed, lost = u_subs(rows[donor])
        candidates, n_contact_diff = sw.null_candidates(rows, pool, donor, columns)
    else:
        donor_sub, placed, lost, placement = {}, 0, 0, {}
        for c in u_idx:
            s_, p_, l_, placement[c] = sw.protein_subs(chains[c]["seq"], ortholog.protein)
            donor_sub.update({(c, i): aa for i, aa in s_.items()}); placed += p_; lost += l_
        ids = cols if cols is not None else {c: list(chains[c]["num"]) for c in u_idx}
        n_contact_diff = sw.protein_contact_differences({c: chains[c]["seq"] for c in u_idx}, placement, ids, flags)
        candidates = sw.null_candidates(rows, pool, donor, columns, n=n_contact_diff, host=host)[0] if rows is not None else []
        if cols is None:
            columns = sorted({n for c in u_idx for n, hit in zip(chains[c]["num"], flags[c]) if hit})
        extra = {**extra, "reading": "protein", "donor_protein": ortholog.accession, "donor_gene_id": ortholog.gene_id}
    al = cybrid.aligner()
    donor_on = {g: cybrid.on_human(mt[g], dp[g][0], dp[g][1], al) for g in genes}
    sets, unmodelled = cs.plan(chains, mt_chain, mt, donor_on, {g: {} for g in genes})
    unmodelled += sum(insertions(mt[g], dp[g][0], al) for g in genes if not dp[g][1])
    mt_sub = sets["all"]
    if not donor_sub:
        null = "not run: the donor's ortholog changes no placed residue, so its restoration is 0 and cannot be a candidate's"
    elif rows is None:
        null = "not assessed: no alignment of the gene (the null draws from the alignments)"
    elif len(candidates) < sw.MIN_CANDIDATES:
        null = f"not assessed: {len(candidates)} candidates, at least {sw.MIN_CANDIDATES} needed"
    else:
        null = "assessed"
    # the pairs of the gene: (mitochondrial chain, chain of u) that touch, whose structure states are Foldseek's for both chains
    states = chain_state_table(chains, structure_tokens(cif, binary, tokens_dir))
    pairs, skipped = [], {"too long": 0, "no structure states": 0}
    for c in u_idx:
        for k, pos in pair_positions(chains, one_idx, c).items():
            if len(chains[k]["seq"]) + LINKER + len(chains[c]["seq"]) > MAX_RESIDUES:
                skipped["too long"] += 1
            elif states[k] is None or states[c] is None:
                skipped["no structure states"] += 1
            else:
                pairs.append((k, c, pos))
    case = sw.case_row(donor, complex_, gene, "scored", host=host, copies=len(u_idx), contact_residues=len(columns), placed=placed, not_placed=lost,
                       substituted=len(donor_sub), contact_differences=n_contact_diff, candidates=len(candidates), null=null,
                       mt_substituted=len(mt_sub), unmodelled=unmodelled, pdb=pid, pairs=len(pairs), pairs_too_long=skipped["too long"],
                       pairs_no_states=skipped["no structure states"], **extra)
    if not pairs:
        why = " and ".join(f"{n} pair(s) with {k}" for k, n in skipped.items() if n) or "no pair of chains touches"
        return done({**case, "status": f"not assessed: no pair of a mitochondrial chain and a chain of the gene can be read ({why})", "null": ""})
    drawn = [sw.draw_species(candidates, k) for k in range(draws)] if null == "assessed" else []
    species_sub = {s: u_subs(rows[s])[0] for s in sorted(set(drawn))}
    if logits_for is None:
        logits_at, _ = load_model(weights_dir, device)
        vocab = read_vocab(weights_dir / "vocab.txt")
    else:
        logits_at, vocab = logits_for(), vocab_for()
    mt_by_chain = {k: {i: aa for (kk, i), aa in mt_sub.items() if kk == k} for k, _, _ in pairs}

    def read(k: int, c: int, pos: list[int], sub: dict[int, str]) -> np.ndarray:
        toks, off = joined(chains[k]["seq"], states[k], chains[c]["seq"], states[c], mt_by_chain[k], sub)
        return masked_logp(logits_at, vocab, toks, pos, batch)

    # every distinct set of u's substitutions is read once, however often it is drawn
    wanted: dict[tuple, dict] = {}
    host_u = {c: chains[c]["seq"] for c in u_idx}
    # a donor whose ortholog changes nothing has B1 = B0: one read, no shuffles, no single changes
    shuf = {j: {(c, i): aa for c in u_idx for i, aa in shuffled(host_u[c], j).items()} for j in range(shuffles)} if donor_sub else {}
    single = {(c, i): {(c, i): aa} for (c, i), aa in sorted(donor_sub.items())}
    for sub in [{}, donor_sub, *species_sub.values(), *shuf.values(), *single.values()]:
        wanted.setdefault(tuple(sorted(sub.items())), sub)
    logp: dict[tuple, dict[tuple[int, int], np.ndarray]] = {}
    for wk, sub in wanted.items():
        logp[wk] = {(k, c): read(k, c, pos, {i: aa for (cc, i), aa in sub.items() if cc == c}) for k, c, pos in pairs}
        print(f"{donor} {complex_} {gene}: {len(logp)} of {len(wanted)} reads of u, {sum(len(p) for _, _, p in pairs)} places each", flush=True)
    key = lambda sub: tuple(sorted(sub.items()))
    total = lambda sub: float(sum(np.nansum(v) for v in logp[key(sub)].values()))
    n_pos = int(sum(np.isfinite(v).sum() for v in logp[key({})].values()))

    def row(kind: str, sub: dict, draw: int, sp: str) -> dict:
        pll = total(sub)
        return {**{k: case[k] for k in ("host", "donor", "complex", "gene")}, "kind": kind, "draw": draw, "species": sp, "seed": SEED,
                "substituted": len(sub), "pairs": len(pairs), "positions": n_pos, "pll": pll, "ddG": -pll}
    res = [row("B0", {}, -1, host), row("B1", donor_sub, -1, donor), *(row("null", species_sub[s], k, s) for k, s in enumerate(drawn)),
           *(row("shuffled", shuf[j], j, "shuffled") for j in sorted(shuf))]
    mt_gene = {k: g for k, (g, _) in mt_chain.items()}
    prow = lambda kind, draw, k, c, v: {**{x: case[x] for x in ("host", "donor", "complex", "gene")}, "kind": kind, "draw": draw, "mt_gene": mt_gene[k],
                                       "mt_chain": chains[k]["id"], "u_chain": chains[c]["id"], "positions": int(np.isfinite(v).sum()), "pll": float(np.nansum(v))}
    prs = [prow(kind, draw, k, c, logp[key(sub)][(k, c)]) for kind, draw, sub in (("B0", -1, {}), ("B1", -1, donor_sub), *(("shuffled", j, shuf[j]) for j in sorted(shuf)))
           for k, c, _ in pairs]
    # the residues the ortholog changes, each alone
    dist = {c: residue_distances(chains, one_idx, c) for c in u_idx}
    resid = []
    for (c, i), aa in sorted(donor_sub.items()):
        for k, cc, pos in pairs:
            if cc != c:
                continue
            delta = logp[key(single[(c, i)])][(k, c)] - logp[key({})][(k, c)]
            fin = np.isfinite(delta)
            top = int(np.nanargmax(np.abs(np.where(fin, delta, 0.0)))) if fin.any() else None
            resid.append({**{x: case[x] for x in ("host", "donor", "complex", "gene")}, "mt_gene": mt_gene[k], "mt_chain": chains[k]["id"],
                          "u_chain": chains[c]["id"], "resnum": chains[c]["num"][i], "host_aa": chains[c]["seq"][i], "donor_aa": aa,
                          "interface": bool(dist[c][i] <= cs.CONTACT_A), "min_dist": float(dist[c][i]), "positions": int(fin.sum()),
                          "delta_pll": float(np.nansum(delta)), "top_mt_resnum": chains[k]["num"][pos[top]] if top is not None else np.nan,
                          "top_mt_delta": float(delta[top]) if top is not None else np.nan})
    cs.report_gpu(device)
    return done(case, pd.DataFrame(res, columns=SCORE_COLUMNS), pd.DataFrame(prs, columns=PAIR_COLUMNS), pd.DataFrame(resid, columns=RESIDUE_COLUMNS))


def probe(gb: Path, pdb_dir: Path, tokens_dir: Path, weights_dir: Path, device: str, pdb: str = "5Z62", mt_chain: str = "A", u_chain: str = "D") -> dict:
    """Read one real pair, the human cytochrome c oxidase's COX1 and COX4I1 (5Z62 chains A and D), as a task reads it, and print what it shows:
    the card, the attention the model runs with, the joined input's length, the vocabulary against the model's own tokenizer, the time and
    memory of a read, the PLL of the contact residues with the real partner and with the partner shuffled, and the same twice. Writes nothing."""
    import time
    import torch
    from transformers import EsmTokenizer
    t0 = time.time()
    logits_at, rev = load_model(weights_dir, device)
    print("model loaded in", round(time.time() - t0, 1), "s;", "cuda" if device.startswith("cuda") else device, torch.cuda.get_device_name(0) if device.startswith("cuda") else "")
    vocab = read_vocab(weights_dir / "vocab.txt")
    cif = pdb_dir / f"{pdb}.cif"
    chains = cs.read_chains(cif)
    tok = structure_tokens(cif, None, tokens_dir)
    states = chain_state_table(chains, tok)
    ix = {c["id"]: i for i, c in enumerate(chains)}
    k, c = ix[mt_chain], ix[u_chain]
    pos = pair_positions(chains, [k], c)[k]
    toks, off = joined(chains[k]["seq"], states[k], chains[c]["seq"], states[c])
    tk = EsmTokenizer.from_pretrained(str(weights_dir))
    text = "".join(toks[:40])
    same = tk(text)["input_ids"] == [vocab["<cls>"], *(vocab[t] for t in toks[:40]), vocab["<eos>"]]
    print(f"joined input {len(toks)} tokens (chain {mt_chain} {len(chains[k]['seq'])}, linker {LINKER}, chain {u_chain} {len(chains[c]['seq'])}), "
          f"{len(pos)} contact places; the first tokens {' '.join(toks[:6])}; the tokenizer gives the same ids as the vocabulary: {same}")
    out = {}
    for name, sub in (("real", {}), ("shuffled", shuffled(chains[c]["seq"], 0))):
        toks, _ = joined(chains[k]["seq"], states[k], chains[c]["seq"], states[c], None, sub)
        t0 = time.time()
        a = masked_logp(logits_at, vocab, toks, pos)
        dt = time.time() - t0
        b = masked_logp(logits_at, vocab, toks, pos)
        out[name] = float(np.nansum(a))
        print(f"{name}: PLL of the {len(pos)} contact places {out[name]:.4f}; {dt:.2f} s a read; the same read again differs by {np.abs(a - b).max():.2e}", flush=True)
    toks, _ = joined(chains[k]["seq"], states[k], chains[c]["seq"], states[c])
    b1 = masked_logp(logits_at, vocab, toks, pos, batch=1)
    print(f"batch of 1 against batch of {BATCH}: largest difference {np.abs(b1 - masked_logp(logits_at, vocab, toks, pos)).max():.2e}")
    cs.report_gpu(device)
    return out


# --- the tables and the rule ----------------------------------------------------------------------------------------


def read_tasks(out: Path, want: list[tuple[str, str, str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    """(cases, scores, pairs, residues, the tasks with no table) from the tasks' tables. "null" is a kind of row, not a missing value."""
    parts: dict[str, list[pd.DataFrame]] = {"cases": [], "": [], "pairs": [], "residues": []}
    missing = []
    for d, c, g in want:
        f = out / STATE / f"{d}.{c}.{g}"
        if not Path(f"{f}.cases.tsv").is_file():
            missing.append(f"{d} / {c} / {g}")
            continue
        for name in parts:
            parts[name].append(pd.read_csv(f"{f}.{name}.tsv" if name else f"{f}.tsv", sep="\t", keep_default_na=False, na_values=[""]))
    if missing:
        print(f"MISSING: {len(missing)} of {len(want)} tasks have no table: " + "; ".join(missing))
    if not parts["cases"]:
        raise SystemExit("no table to collect")
    cat = {k: pd.concat(v, ignore_index=True) for k, v in parts.items()}
    # a task with no row leaves an empty table whose columns read as text, and the join then holds numbers and booleans as objects
    for name, cols in (("", NUMERIC_SCORES), ("pairs", NUMERIC_PAIRS), ("residues", NUMERIC_RESIDUES)):
        for c in cols:
            cat[name][c] = pd.to_numeric(cat[name][c])
    cat["residues"]["interface"] = cat["residues"]["interface"].astype(str).eq("True")
    return cat["cases"], cat[""], cat["pairs"], cat["residues"], missing


def dependence_table(res: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """One row per scored (donor, complex, gene): PLL(B0) minus the mean PLL of the shuffled partners (positive: the real partner makes the
    mitochondrial residues likelier than a shuffled one; 0 would be a model that does not use the partner) and the number of shuffles."""
    rows = []
    for r in res.itertuples():
        s = scores[(scores.donor == r.donor) & (scores["complex"] == r.complex) & (scores.gene == r.gene)]
        sh = s[s.kind == "shuffled"]
        b0 = s[s.kind == "B0"].pll.iloc[0]
        rows.append({"donor": r.donor, "complex": r.complex, "gene": r.gene, "pll_b0": b0, "pll_shuffled_mean": sh.pll.mean() if len(sh) else np.nan,
                     "real_minus_shuffled": b0 - sh.pll.mean() if len(sh) else np.nan, "shuffles": len(sh)})
    return pd.DataFrame(rows, columns=["donor", "complex", "gene", "pll_b0", "pll_shuffled_mean", "real_minus_shuffled", "shuffles"])


def residue_table(residues: pd.DataFrame) -> pd.DataFrame:
    """One row per changed residue of a nuclear chain in a (donor, complex, gene): the change of the PLL summed over the mitochondrial
    chains it faces, the residue's distance to them and whether it is at the interface; and its rank among the gene's interface residues."""
    if not len(residues):
        return pd.DataFrame(columns=["donor", "complex", "gene", "u_chain", "resnum", "host_aa", "donor_aa", "interface", "min_dist", "delta_pll", "rank_in_gene"])
    key = ["donor", "complex", "gene", "u_chain", "resnum", "host_aa", "donor_aa"]
    g = residues.groupby(key, sort=False)
    t = pd.concat([g.interface.first(), g.min_dist.first(), g.delta_pll.sum()], axis=1).reset_index()
    i = t.interface.astype(bool)
    t["rank_in_gene"] = np.nan
    t.loc[i, "rank_in_gene"] = t[i].groupby(["donor", "complex", "gene"]).delta_pll.rank(method="min", ascending=False)
    return t


def interface_contrast(t: pd.DataFrame) -> dict:
    """Does the model's reading of a changed residue depend on the residue touching the mitochondrial chains? The absolute change of the PLL at
    changed interface residues against changed residues away from them (Mann-Whitney, one-sided: the interface larger), and the
    Spearman correlation of the absolute change with the distance."""
    a = t.delta_pll.abs()
    if not len(t) or t.interface.nunique() < 2:
        return {"n_interface": int(t.interface.sum()) if len(t) else 0, "n_away": int((~t.interface).sum()) if len(t) else 0, "u_p": np.nan, "auc": np.nan, "rho_dist": np.nan, "p_dist": np.nan}
    x, y = a[t.interface], a[~t.interface]
    u = stats.mannwhitneyu(x, y, alternative="greater")
    rho, p = stats.spearmanr(t.min_dist, a)
    return {"n_interface": len(x), "n_away": len(y), "u_p": u.pvalue, "auc": u.statistic / (len(x) * len(y)), "rho_dist": rho, "p_dist": p}


def collect(out: Path, measurements: Path = cybrid.MEASUREMENTS, host: str = sw.HOST) -> pd.DataFrame:
    """The tables of the step-2 rule on SaProt's reads (saprot_scores, _pairs, _cases, _restoration, _candidates, _sums, _spearman, and for
    another host _unassessed and _reduced), the per-residue tables (saprot_residues, saprot_residue_summary) and the partner's use
    (saprot_dependence). A task without a table is listed as missing, never as zero. The rule's code is the swap's, applied to the
    column `ddG`, which holds the misfit (minus the PLL): B0 minus B1 is then PLL(B1) minus PLL(B0)."""
    want = sw.tasks(measurements, host)
    cases, scores, pairs, residues, _ = read_tasks(out, want)
    res = sw.restoration_table(cases, scores)
    cand = sw.candidate_table(cases, res)
    meas = cybrid.read_measurements(measurements)
    sums, rho = sw.activity_table(res, cases, meas, host)
    rt = residue_table(residues)
    tables = [("scores", scores), ("pairs", pairs), ("cases", cases), ("restoration", res), ("candidates", cand), ("sums", sums), ("spearman", rho),
              ("dependence", dependence_table(res, scores)), ("residues", rt), ("residue_contrast", pd.DataFrame([interface_contrast(rt)]))]
    if host != sw.HOST:
        tables += [("unassessed", sw.unassessed_donors(measurements, host)), ("reduced", sw.reduced_table(res, cases, meas, host))]
    for name, t in tables:
        t.to_csv(out / f"saprot_{name}.tsv", sep="\t", index=False, float_format="%.5f")
    with pd.option_context("display.width", 250, "display.max_rows", 400, "display.max_columns", 40):
        print(res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(cand[cand.donors_with_ortholog > 0].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(sums.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(rho.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    return res


# --- the bar: the comparison with the count (step 1) and StaB-ddG (step 2) on the same cases -----------------------------


def stratified_auc(df: pd.DataFrame, col: str) -> tuple[float, int]:
    """(the share of (reduced, not reduced) pairs of cases of one host and complex in which `col` is larger in the reduced one, ties
    counting half; the number of such pairs). Cases of different complexes or hosts are never paired."""
    wins, n = 0.0, 0
    for _, g in df.groupby(["host", "complex"]):
        r, nr = g[g.status == "reduced"][col].to_numpy(), g[g.status == "not_reduced"][col].to_numpy()
        for x in r:
            wins += float(((x > nr) + 0.5 * (x == nr)).sum()); n += len(nr)
    return (wins / n if n else np.nan), n


def compare(saprot: dict[str, Path], stab: dict[str, Path], exposure: Path, measurements: Path = cybrid.MEASUREMENTS) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(one row per (host, donor, complex) with a stated outcome, the metrics): the three readouts summed over the same genes, SaProt's
    restoration, StaB-ddG's (step 2) and the count of the residues of those subunits that face a donor difference (step 1), each larger for
    a cybrid predicted worse. The genes are those that SaProt and StaB-ddG both scored for every donor of the complex of the host. The
    metrics: the stratified AUC and the Spearman correlation with the stated activity (expected negative), per readout.
    `saprot[host]` and `stab[host]` are the folders of the tables of that host; `exposure` is step 1's nuclear_exposure.tsv."""
    from .nuclear_where import labelled
    expo = pd.read_csv(exposure, sep="\t")
    meas = cybrid.read_measurements(measurements)
    rows = []
    for host, sdir in saprot.items():
        sap = pd.read_csv(sdir / "saprot_restoration.tsv", sep="\t")
        stb = pd.read_csv(stab[host] / "swap_restoration.tsv", sep="\t")
        lab = labelled(meas[meas.host == host], complexes=tuple(COMPLEXES))
        for cx, g in sap.groupby("complex"):
            donors = sorted(set(g.donor))
            common = set.intersection(*(set(g[g.donor == d].gene) & set(stb[(stb.donor == d) & (stb["complex"] == cx)].gene) for d in donors))
            for d in donors:
                e = expo[(expo.host == host) & (expo.donor == d) & (expo["complex"] == cx)]
                have = [x for x in common if x in e.columns]
                rows.append({"host": host, "donor": d, "complex": cx, "genes": len(common),
                             "saprot": g[(g.donor == d) & g.gene.isin(common)].restoration.sum(),
                             "stabddg": stb[(stb.donor == d) & (stb["complex"] == cx) & stb.gene.isin(common)].restoration.sum(),
                             "count": float(e[have].sum(axis=1).iloc[0]) if len(e) and have else np.nan})
    cases = pd.DataFrame(rows).merge(pd.concat([labelled(meas[meas.host == h], complexes=tuple(COMPLEXES)) for h in saprot])[["host", "donor", "complex", "status", "pct"]],
                                      on=["host", "donor", "complex"], how="left")
    lab = cases[cases.status.isin(["reduced", "not_reduced"])]
    out = []
    for name in ("saprot", "stabddg", "count"):
        auc, n = stratified_auc(lab.dropna(subset=[name]), name)
        s = lab.dropna(subset=["pct", name])
        rho, p = stats.spearmanr(s.pct, s[name]) if len(s) >= 3 and s.pct.nunique() > 1 and s[name].nunique() > 1 else (np.nan, np.nan)
        out.append({"readout": name, "cases": len(lab.dropna(subset=[name])), "pairs": n, "stratified_auc": auc, "n_activity": len(s), "spearman_activity": rho, "p": p})
    return cases, pd.DataFrame(out)
