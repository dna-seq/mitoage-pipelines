"""Downloaders for everything the analyses read. Each download appends a record to data/input/MANIFEST.json.

The MitoAge build: the Hugging Face dataset `longevity-genie/mitoage` (tables and alignments).

NCBI ortholog sets: what `datasets download gene symbol <GENE> --taxon human --ortholog mammalia
--include protein,cds` gives, through the NCBI Datasets v2 REST API in three calls:
  GET  /gene/symbol/<GENE>/taxon/human                    the human gene's id
  GET  /gene/id/<id>/orthologs?taxon_filter=mammalia      the gene ids of its ortholog set (paged; the
                                                          human gene itself is among them)
  POST /gene/download                                     the zip: protein.faa, cds.fna, data_report.jsonl
The proxy variables of the environment (https_proxy) are honoured. NCBI_API_KEY, when set, is sent.
"""
import io
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

HF_REPO = "longevity-genie/mitoage"
NCBI_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"


# --- manifest ---------------------------------------------------------------------------------
def folder_stats(folder: Path) -> tuple[int, int]:
    """(number of files, total bytes) under `folder`, the downloader's own cache folders left out."""
    files = [p for p in sorted(folder.rglob("*")) if p.is_file() and ".cache" not in p.relative_to(folder).parts]
    return len(files), sum(p.stat().st_size for p in files)


def append_manifest(manifest: Path, source: str, identifier: str, folder: Path, **extra) -> dict:
    """Add one record (never rewrite an earlier one): source, identifier, UTC time, file count, bytes."""
    records = json.loads(manifest.read_text()) if manifest.is_file() else []
    if not isinstance(records, list):
        raise ValueError(f"{manifest} is not a list of download records")
    n, size = folder_stats(folder)
    record = {"source": source, "identifier": identifier,
              "downloaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
              "files": n, "bytes": size, **extra}
    records.append(record)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(records, indent=2) + "\n")
    return record


# --- prepared data ----------------------------------------------------------------------------
# What we fetched or computed, published by topic, each dataset laid out like data/ (input/, interim/, output/).
PREPARED = {
    "partners": "longevity-genie/mitoage-nuclear-partners",    # NCBI ortholog sets as fetched (mammals; complex IV across vertebrates); mitoport's alignments and trees
    "embeddings": "longevity-genie/mitoage-embeddings",        # ESM-2 650M mean embeddings of both partners
    "residues": "longevity-genie/mitoage-embeddings",          # the same, per residue (7.7 GB): what `mitoport iface` reads
    "structures": "longevity-genie/mitoage-structures",        # 806 ESMFold2 structures of the 13 proteins in 62 mammals
    "aligners": "longevity-genie/mitoage-aligner-variants",    # the 13 proteins aligned by four aligners
    "results": "longevity-genie/mitoage-pipelines",            # the dated tree and the result tables, complex IV's contact map among them
}


def prepared(topic: str, root: Path, manifest: Path, revision: str | None = None, results_dir: Path | None = None) -> None:
    """One topical dataset into the data root. NCBI's zips are unpacked where `download orthologs` puts them.
    The `results` topic holds our result tables: they go to `results_dir` (data/reference), never over
    data/output, so that a rerun can be compared with them."""
    from huggingface_hub import HfApi, snapshot_download
    repo = PREPARED[topic]
    sha = HfApi().dataset_info(repo, revision=revision).sha
    if topic == "results":
        dest = results_dir or root / "reference"
        snapshot_download(repo, repo_type="dataset", revision=sha, local_dir=dest, allow_patterns=["output/*", "output/**"])
        snapshot_download(repo, repo_type="dataset", revision=sha, local_dir=root, allow_patterns=["input/**", "interim/**"])
    else:
        dest = root
        # the per-residue embeddings share a dataset with the means: only the `residues` topic fetches them
        allow = ["input/embeddings/*.residues.npy", "input/embeddings/*.index.parquet"] if topic == "residues" else ["input/**", "interim/**"]
        ignore = ["*.residues.npy"] if topic == "embeddings" else None
        snapshot_download(repo, repo_type="dataset", revision=sha, local_dir=root, allow_patterns=allow, ignore_patterns=ignore)
    if topic == "partners":           # the mammalian sets, and complex IV's eleven subunits across vertebrates (H12)
        for z in sorted([*(root / "input" / "ncbi").glob("*.zip"), *(root / "input" / "ncbi-vertebrata").glob("*.zip")]):
            gene_dir = z.with_suffix("")
            if not (gene_dir / "ncbi_dataset" / "data" / "protein.faa").is_file():
                gene_dir.mkdir(exist_ok=True)
                with zipfile.ZipFile(z) as f: f.extractall(gene_dir)
    record = append_manifest(manifest, "huggingface", f"datasets/{repo}", dest, revision=sha, topic=topic)
    print(f"{repo} at {sha} ({topic}) into {dest}")


# --- the dated species tree -------------------------------------------------------------------
SITE_REPO = "dna-seq/mitoage"
SITE_TREE_PATH = "src/data/generated/phylogeny.json"


def site_tree(dest: Path, manifest: Path, ref: str = "main") -> Path:
    """The build's dated species tree (Open Tree topology, TimeTree dates; a leaf is a species id) from the
    site's repository, at the commit `ref` points to. It is not part of the Hugging Face dataset."""
    sha = json.loads(_request(f"https://api.github.com/repos/{SITE_REPO}/commits/{ref}"))["sha"]
    body = _request(f"https://raw.githubusercontent.com/{SITE_REPO}/{sha}/{SITE_TREE_PATH}", accept="*/*")
    json.loads(body)                                    # refuse anything that is not the tree
    dest.parent.mkdir(parents=True, exist_ok=True); dest.write_bytes(body)
    record = append_manifest(manifest, "github", f"{SITE_REPO}/{SITE_TREE_PATH}", dest.parent, revision=sha)
    print(f"{SITE_TREE_PATH} of {SITE_REPO} at {sha}: {record['bytes']} bytes in {dest}")
    return dest


# --- the MitoAge build ------------------------------------------------------------------------
def build(dest: Path, manifest: Path, revision: str | None = None) -> Path:
    """The Hugging Face dataset into `dest`. `revision`: a commit of the dataset, to pin what is read."""
    from huggingface_hub import HfApi, snapshot_download
    sha = HfApi().dataset_info(HF_REPO, revision=revision).sha      # the commit actually downloaded
    dest.mkdir(parents=True, exist_ok=True)
    snapshot_download(HF_REPO, repo_type="dataset", revision=sha, local_dir=dest)
    record = append_manifest(manifest, "huggingface", f"datasets/{HF_REPO}", dest, revision=sha)
    print(f"{HF_REPO} at {sha}: {record['files']} files, {record['bytes']} bytes in {dest}")
    return dest


# --- solved structures -----------------------------------------------------------------------
RCSB = "https://files.rcsb.org/download"


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def genbank(accessions: list[str], dest: Path, manifest: Path) -> Path:
    """GenBank records (with their CDS translations) by accession.version, through E-utilities efetch, into
    `dest`/records.gb, with the file's SHA-256 in the manifest. A versioned accession names fixed bytes of
    sequence; the annotation around it can still change, which the digest records."""
    import hashlib
    dest.mkdir(parents=True, exist_ok=True)
    unversioned = [a for a in accessions if "." not in a]
    if unversioned:
        raise ValueError(f"accessions without a version: {unversioned}")
    ids = ",".join(sorted(set(accessions)))
    out = dest / "records.gb"
    out.write_bytes(_request(f"{EUTILS}/efetch.fcgi?db=nuccore&id={ids}&rettype=gbwithparts&retmode=text", accept="*/*"))
    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    append_manifest(manifest, "ncbi-efetch", ids, dest, sha256=digest)
    print(f"{out}  sha256 {digest}")
    return out


def pdb(ids: list[str], dest: Path, manifest: Path) -> list[Path]:
    """mmCIF files of solved structures from RCSB into `dest`, one record per call with each file's SHA-256:
    RCSB revises entries in place, so the digest is what says which file a result was computed on.
    A file already there is kept and not fetched again."""
    import hashlib
    dest.mkdir(parents=True, exist_ok=True); got, digests = [], {}
    for i in ids:
        p = dest / f"{i}.cif"
        if not p.is_file():
            p.write_bytes(_request(f"{RCSB}/{i}.cif", accept="*/*"))
        digests[i] = hashlib.sha256(p.read_bytes()).hexdigest(); got.append(p)
    append_manifest(manifest, "rcsb", ",".join(ids), dest, sha256=digests)
    for i, d in digests.items(): print(f"{i}.cif  sha256 {d}")
    return got


def stabddg(dest: Path, manifest: Path) -> list[Path]:
    """StaB-ddG's three weight files (H28) from its repository at the pinned commit into `dest`. Each is kept only
    when its SHA-256 is the pinned one; a file already there with that digest is not fetched again."""
    import hashlib
    from .mitoport.cybrid_stab import COMMIT, REPO, WEIGHTS
    dest.mkdir(parents=True, exist_ok=True); got = []
    raw = REPO.replace("https://github.com/", "https://raw.githubusercontent.com/")
    for name, (file, digest) in WEIGHTS.items():
        p = dest / file
        if not (p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest() == digest):
            body = _request(f"{raw}/{COMMIT}/model_ckpts/{file}", accept="*/*")
            if hashlib.sha256(body).hexdigest() != digest:
                raise SystemExit(f"{file} at {COMMIT[:7]} is not the pinned file (SHA-256 differs): nothing written")
            p.write_bytes(body)
        got.append(p); print(f"{file}  sha256 {digest}")
    append_manifest(manifest, "github", f"{REPO}@{COMMIT}", dest, sha256={f: d for f, d in WEIGHTS.values()})
    return got


# --- NCBI ortholog sets -----------------------------------------------------------------------
def _request(url: str, body: dict | None = None, accept: str = "application/json", tries: int = 5) -> bytes:
    headers = {"accept": accept, "user-agent": "mitoage-pipelines"}
    if os.environ.get("NCBI_API_KEY"): headers["api-key"] = os.environ["NCBI_API_KEY"]
    data = None
    if body is not None:
        data = json.dumps(body).encode(); headers["content-type"] = "application/json"
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=headers), timeout=300) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or attempt == tries - 1:
                raise RuntimeError(f"the server answered {e.code} for {url}: {e.read()[:300]!r}") from e
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == tries - 1:
                raise RuntimeError(f"no answer from {url}: {e} (behind a proxy, set https_proxy)") from e
        time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


# genes NCBI has renamed since the analysis named them: our name (folders, tables) -> NCBI's current symbol.
# NDUFA4 is COXFA4 at NCBI (gene 4697) as of 2026-11-15; the old symbol now matches no gene exactly.
NCBI_SYMBOL = {"NDUFA4": "COXFA4"}


def ortholog_gene_ids(symbol: str, taxon: str = "human", ortholog_taxon: str = "mammalia") -> tuple[str, list[int]]:
    """(id of the gene in `taxon`, ids of its ortholog set within `ortholog_taxon`, in NCBI's order). `symbol` is our
    name; a gene NCBI has renamed is looked up by its current symbol (NCBI_SYMBOL), and the substitution printed."""
    q = urllib.parse.quote
    ncbi = NCBI_SYMBOL.get(symbol, symbol)
    if ncbi != symbol:
        print(f"{symbol}: looked up at NCBI as {ncbi}, its current symbol")
    found = json.loads(_request(f"{NCBI_API}/gene/symbol/{q(ncbi)}/taxon/{q(taxon)}")).get("reports", [])
    ids = [r["gene"]["gene_id"] for r in found if r["gene"]["symbol"] == ncbi]
    if len(ids) != 1:
        raise RuntimeError(f"{symbol} in {taxon}: expected one gene at NCBI, found {len(ids)}")
    gene_ids: list[int] = []; token = None
    while True:
        url = f"{NCBI_API}/gene/id/{ids[0]}/orthologs?taxon_filter={q(ortholog_taxon)}&page_size=1000"
        if token: url += f"&page_token={q(token)}"
        page = json.loads(_request(url))
        gene_ids += [int(r["gene"]["gene_id"]) for r in page.get("reports", [])]
        token = page.get("next_page_token")
        if not token: break
    if len(set(gene_ids)) != len(gene_ids):
        raise RuntimeError(f"{symbol}: NCBI listed a gene id twice in the ortholog set")
    if not gene_ids:
        raise RuntimeError(f"{symbol}: NCBI lists no orthologs in {ortholog_taxon}")
    return ids[0], gene_ids


def orthologs_of(symbol: str, dest: Path, manifest: Path, taxon: str = "human", ortholog_taxon: str = "mammalia") -> Path:
    """The ortholog set of one gene, unzipped into `dest`/<GENE>/ (ncbi_dataset/data/protein.faa, cds.fna, ...)."""
    human_id, gene_ids = ortholog_gene_ids(symbol, taxon, ortholog_taxon)
    blob = _request(f"{NCBI_API}/gene/download", accept="application/zip",
                    body={"gene_ids": gene_ids, "include_annotation_type": ["FASTA_PROTEIN", "FASTA_CDS"]})
    folder = dest / symbol
    folder.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        for name in z.namelist():           # refuse a member that would land outside the folder
            if not (folder / name).resolve().is_relative_to(folder.resolve()):
                raise RuntimeError(f"{symbol}: the archive holds a path outside its folder: {name}")
        z.extractall(folder)
    data = folder / "ncbi_dataset" / "data"
    for needed in ("protein.faa", "cds.fna"):
        if not (data / needed).is_file():
            raise RuntimeError(f"{symbol}: the download holds no {needed}")
    record = append_manifest(
        manifest, "ncbi-datasets-v2", f"gene/symbol/{symbol}/taxon/{taxon}/orthologs/{ortholog_taxon}", folder,
        gene_id=human_id, ortholog_genes=len(gene_ids), include=["protein", "cds"])
    print(f"{symbol}: {len(gene_ids)} genes in the ortholog set, {record['files']} files, {record['bytes']} bytes")
    return folder


def orthologs(symbols: list[str], dest: Path, manifest: Path, taxon: str = "human",
              ortholog_taxon: str = "mammalia", skip_existing: bool = True) -> None:
    """The ortholog sets of every gene in `symbols`. A gene already downloaded is skipped, and said so."""
    skipped = []
    for symbol in symbols:
        if skip_existing and (dest / symbol / "ncbi_dataset" / "data" / "protein.faa").is_file():
            skipped.append(symbol); continue
        orthologs_of(symbol, dest, manifest, taxon, ortholog_taxon)
    if skipped:
        print(f"already there, not downloaded again ({len(skipped)}): {' '.join(skipped)}")
