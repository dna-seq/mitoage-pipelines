"""Downloaders for every source of build 2. Each appends one record to data/input/MANIFEST.json (source,
identifier, revision or version, UTC time, files, bytes, md5 where the source has no version) and refuses,
with one line, a file that does not match what it should be.

NCBI E-utilities: at most 3 requests a second, batches of 100 versions, NCBI_API_KEY sent when set, and a
contact address only when NCBI_EMAIL is set (none is written in the code).
"""
import csv
import gzip
import hashlib
import html
import io
import json
import os
import re
import shutil
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

from .. import download as base
from . import pins

UA = "mitoage-pipelines"


class Refused(RuntimeError):
    """A download that is not what it should be: wrong checksum, wrong release, missing records."""


def write_label(folder: Path, label: str) -> None:
    """How the steps name this input in the tables (changelog upstream_version, reasons)."""
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "label.json").write_text(json.dumps({"label": label}) + "\n")


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _open(url: str, data: bytes | None = None, headers: dict | None = None, tries: int = 5, method: str | None = None):
    hdr = {"user-agent": UA, **(headers or {})}
    for attempt in range(tries):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, data=data, headers=hdr, method=method), timeout=600)
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or attempt == tries - 1:
                raise RuntimeError(f"HTTP {e.code} for {url}") from e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            if attempt == tries - 1:
                raise RuntimeError(f"no answer from {url}: {e} (behind a proxy, set https_proxy)") from e
        time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def fetch(url: str, **kw) -> bytes:
    with _open(url, **kw) as r:
        return r.read()


def fetch_to(url: str, dest: Path) -> dict:
    """Stream `url` into `dest`; returns the response headers that date it."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    with _open(url) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 20)
        headers = {k: r.headers.get(k) for k in ("Last-Modified", "Content-Length", "ETag") if r.headers.get(k)}
    os.replace(tmp, dest)
    return headers


def _manifest(manifest: Path, source: str, identifier: str, folder: Path, **extra) -> dict:
    rec = base.append_manifest(manifest, source, identifier, folder, **extra)
    print(f"{source} {identifier}: {rec['files']} files, {rec['bytes']} bytes in {folder}")
    return rec


# --- Hugging Face: build 1 and the published build 2 -------------------------------------------
def hf(dest: Path, manifest: Path, revision: str = pins.HF_BUILD2_REVISION) -> Path:
    """build1/, the build 2 tables and the alignment FASTA files of longevity-genie/mitoage at `revision`."""
    from huggingface_hub import HfApi, snapshot_download
    sha = HfApi().dataset_info(pins.HF_REPO, revision=revision).sha
    if not sha.startswith(revision):
        raise Refused(f"{pins.HF_REPO}: revision {revision} resolved to {sha}")
    dest.mkdir(parents=True, exist_ok=True)
    snapshot_download(pins.HF_REPO, repo_type="dataset", revision=sha, local_dir=dest,
                      allow_patterns=["build1/*", "*.parquet", "review_list.csv", "README.md",
                                      "alignments/*.fasta.gz", "alignments/alignments.json"],
                      ignore_patterns=["alignments/alignments.parquet"])
    if not (dest / "build1" / "total_mtDNA_base_composition.parquet").is_file():
        raise Refused(f"{pins.HF_REPO} at {sha} has no build1/total_mtDNA_base_composition.parquet")
    _manifest(manifest, "huggingface", f"datasets/{pins.HF_REPO}", dest, revision=sha,
              matches_build2_pin=sha == pins.HF_BUILD2_REVISION)
    return dest


# --- RefSeq mitochondrion release -------------------------------------------------------------
def refseq_current(dest_root: Path, manifest: Path, expect: str | None = None) -> Path:
    """The current release from NCBI's FTP site, with its release number and file date."""
    release = fetch(pins.REFSEQ_RELEASE_URL).decode().strip()
    if expect and release != expect:
        raise Refused(f"NCBI's current RefSeq release is {release}, not {expect}: NCBI keeps only the current "
                      f"release; for release {expect} use `download refseq --from-file` with a kept copy, or --efetch")
    dest = dest_root / f"release{release}" / "mitochondrion.1.genomic.gbff.gz"
    headers = fetch_to(pins.REFSEQ_URL, dest)
    (dest.parent / "RELEASE_NUMBER").write_text(release + "\n")
    from email.utils import parsedate_to_datetime
    fdate = parsedate_to_datetime(headers["Last-Modified"]).date().isoformat() if headers.get("Last-Modified") else "unknown"
    write_label(dest.parent, f"release {release} (mitochondrion.1.genomic.gbff.gz, file date {fdate})")
    md5 = md5_of(dest)
    _manifest(manifest, "ncbi-ftp", pins.REFSEQ_URL, dest.parent, revision=f"release {release}",
              last_modified=headers.get("Last-Modified"), md5=md5,
              matches_build2_pin=md5 == pins.REFSEQ_MD5)
    if md5 != pins.REFSEQ_MD5:
        print(f"NOTE: release {release} is not build 2's release {pins.REFSEQ_RELEASE} (md5 {md5} != {pins.REFSEQ_MD5})")
    return dest


def refseq_from_file(src: Path, dest_root: Path, manifest: Path, release: str = pins.REFSEQ_RELEASE,
                     md5: str = pins.REFSEQ_MD5) -> Path:
    """A kept copy of a past release, accepted only when its md5 is the pinned one."""
    if not src.is_file():
        raise FileNotFoundError(f"no file at {src}")
    got = md5_of(src)
    if got != md5:
        raise Refused(f"{src} has md5 {got}; release {release} of build 2 is {md5}")
    dest = dest_root / f"release{release}" / "mitochondrion.1.genomic.gbff.gz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dest)
    (dest.parent / "RELEASE_NUMBER").write_text(release + "\n")
    write_label(dest.parent, pins.REFSEQ_VERSION)
    _manifest(manifest, "local-copy", str(src), dest.parent, revision=f"release {release}", md5=got,
              origin=pins.REFSEQ_URL, matches_build2_pin=True)
    return dest


def _efetch_batch(ids: list[str]) -> str:
    form = {"db": "nuccore", "rettype": "gbwithparts", "retmode": "text", "tool": UA, "id": ",".join(ids)}
    if os.environ.get("NCBI_API_KEY"): form["api_key"] = os.environ["NCBI_API_KEY"]
    if os.environ.get("NCBI_EMAIL"): form["email"] = os.environ["NCBI_EMAIL"]
    body = urllib.parse.urlencode(form).encode()
    return fetch(pins.EFETCH_URL, data=body, headers={"content-type": "application/x-www-form-urlencoded"}).decode()


def split_genbank(text: str) -> list[str]:
    """The records of a multi-record GenBank text, each ending with `//`."""
    out = []
    for chunk in re.split(r"^//\s*$", text, flags=re.M):
        chunk = chunk.strip()
        if chunk.startswith("LOCUS"):
            out.append(chunk + "\n//\n")
    return out


def version_of(record: str) -> str | None:
    m = re.search(r"^VERSION\s+(\S+)", record, flags=re.M)
    return m.group(1) if m else None


def _have(folder: Path, v: str) -> bool:
    """A record of `v` is in `folder` (`v` with a version, or an accession whose file carries the version)."""
    return (folder / f"{v}.gb").is_file() or ("." not in v and any(folder.glob(f"{v}.*.gb")))


def efetch(versions: list[str], folder: Path, batch: int = 100, per_second: float = 3.0) -> list[str]:
    """GenBank records (gbwithparts) of the given accession.versions (or accessions: then NCBI's current
    version), one file per version in `folder`. Ids already there are not fetched again. Returns the ids
    NCBI did not return."""
    folder.mkdir(parents=True, exist_ok=True)
    want = [v for v in dict.fromkeys(versions) if not _have(folder, v)]
    gap, last = 1.0 / per_second, 0.0
    for i in range(0, len(want), batch):
        ids = want[i:i + batch]
        for attempt in range(3):
            time.sleep(max(0.0, last + gap - time.monotonic())); last = time.monotonic()
            got = {version_of(r): r for r in split_genbank(_efetch_batch(ids))}
            for v in ids:
                hit = v if v in got else next((k for k in got if k and k.split(".")[0] == v and "." not in v), None)
                if hit: (folder / f"{hit}.gb").write_text(got[hit])
            ids = [v for v in ids if not _have(folder, v)]
            if not ids: break
        print(f"  efetch {min(i + batch, len(want))}/{len(want)}", flush=True)
    return [v for v in dict.fromkeys(versions) if not _have(folder, v)]


def genbank(versions: list[str], folder: Path, manifest: Path, what: str) -> Path:
    missing = efetch(versions, folder)
    if missing:
        raise Refused(f"NCBI efetch returned no record for {len(missing)} of {len(versions)} versions: {', '.join(missing[:10])}")
    _manifest(manifest, "ncbi-efetch", f"nuccore gbwithparts: {what}", folder, records=len(set(versions)),
              revision="accession.version (efetch returns NCBI's current text of each version)")
    return folder


def refseq_efetch(versions: list[str], dest_root: Path, manifest: Path, label: str) -> Path:
    """Pinned mode without the release file: exactly the given accession.versions, by efetch, as one gzipped file."""
    folder = dest_root / f"efetch-{label}" / "records"
    missing = efetch(versions, folder)
    if missing:
        raise Refused(f"efetch returned no record for {len(missing)} versions: {', '.join(missing[:10])}")
    out = dest_root / f"efetch-{label}" / "mitochondrion.efetch.gbff.gz"
    write_label(out.parent, f"the {len(set(versions))} RefSeq accession.versions of build 2 by E-utilities efetch (retrieved {time.strftime('%Y-%m-%d', time.gmtime())})")
    with gzip.open(out, "wt") as f:
        for v in dict.fromkeys(versions):
            f.write((folder / f"{v}.gb").read_text())
    _manifest(manifest, "ncbi-efetch", f"nuccore gbwithparts: the {len(set(versions))} RefSeq versions of build 2",
              out.parent, revision=f"version list of {label}", md5=md5_of(out))
    return out


# --- NCBI Taxonomy ----------------------------------------------------------------------------
def _extract(names: list[str], read, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    for n in pins.TAXDUMP_FILES:
        if n not in names:
            raise Refused(f"the taxonomy dump has no {n}")
        with read(n) as src, open(dest / n, "wb") as out:
            shutil.copyfileobj(src, out, 1 << 20)


def taxdump_archive(date: str, dest_root: Path, manifest: Path) -> Path:
    """The monthly archive new_taxdump_<date>.zip, checked against its .md5 file."""
    url = pins.TAXDUMP_ARCHIVE_URL.format(date=date)
    dest = dest_root / date
    z = dest / f"new_taxdump_{date}.zip"
    fetch_to(url, z)
    got = md5_of(z)
    try:                                   # the archive publishes no .md5 file (checked 2026-11-14); use it when there is one
        want = fetch(url + ".md5", tries=1).decode().split()[0]
    except RuntimeError:
        want = None
    if want and got != want:
        raise Refused(f"{url}: md5 {got}, NCBI's .md5 file says {want}")
    with zipfile.ZipFile(z) as zf:
        zf.testzip()                       # the zip's own CRCs are the check when no .md5 is published
        names = {Path(n).name: n for n in zf.namelist()}
        _extract(list(names), lambda n: zf.open(names[n]), dest)
    z.unlink()
    write_label(dest, f"new_taxdump archive of {date}")
    _manifest(manifest, "ncbi-ftp", url, dest, revision=f"new_taxdump archive of {date}", md5=got,
              md5_published=want or "none (the archive has no .md5 files; zip CRCs checked)",
              matches_build2_pin="no (build 2 read the current dump of 2026-11-13; closest archive)" if date != pins.TAXDUMP_BUILD2 else "yes")
    return dest


def taxdump_current(dest_root: Path, manifest: Path) -> Path:
    url = pins.TAXDUMP_CURRENT_URL
    tmp = dest_root / "new_taxdump.tar.gz"
    headers = fetch_to(url, tmp)
    want = fetch(url + ".md5").decode().split()[0]
    got = md5_of(tmp)
    if got != want:
        raise Refused(f"{url}: md5 {got}, NCBI's .md5 file says {want}")
    from email.utils import parsedate_to_datetime
    date = parsedate_to_datetime(headers["Last-Modified"]).date().isoformat()
    dest = dest_root / f"current-{date}"
    with tarfile.open(tmp) as tf:
        _extract(tf.getnames(), lambda n: tf.extractfile(n), dest)
    tmp.unlink()
    lm = parsedate_to_datetime(headers["Last-Modified"]).strftime("%Y-%m-%dT%H:%M:%SZ")
    write_label(dest, f"new_taxdump of {date} (Last-Modified {lm})")
    _manifest(manifest, "ncbi-ftp", url, dest, revision=f"new_taxdump of {date}", last_modified=headers["Last-Modified"],
              md5=got, matches_build2_pin=got == pins.TAXDUMP_BUILD2_MD5)
    return dest


# --- lifespan sources -------------------------------------------------------------------------
def anage(dest: Path, manifest: Path) -> Path:
    fetch_to(pins.ANAGE_URL, dest)
    got = md5_of(dest)
    same = got == pins.ANAGE_MD5
    _manifest(manifest, "hagr", pins.ANAGE_URL, dest.parent, md5=got, matches_build2_pin=same,
              revision="unversioned URL")
    if not same:
        print(f"WARNING: AnAge dataset.zip has md5 {got}, build 2 read {pins.ANAGE_MD5} ({pins.ANAGE_VERSION}). "
              f"AnAge has changed since build 2: lifespans will differ.")
    return dest


def fishbase(db: str, dest_root: Path, manifest: Path, version: str = pins.FISHBASE_VERSION) -> Path:
    """The tables the rules read, from source.coop, each checked against build 2's md5 when the version is pinned."""
    dest = dest_root / version
    md5s = {}
    for t in pins.FISHBASE_TABLES:
        url = pins.SOURCECOOP_URL.format(db=db, version=version, table=t)
        p = dest / f"{t}.parquet"
        fetch_to(url, p)
        md5s[t] = md5_of(p)
        pin = pins.FISHBASE_MD5[db].get(t) if version == pins.FISHBASE_VERSION else None
        if pin and md5s[t] != pin:
            raise Refused(f"{url}: md5 {md5s[t]}, build 2 read {pin}")
    name = "fishbase" if db == "fb" else "sealifebase"
    _manifest(manifest, "source.coop", f"cboettig/fishbase/{db}/{version}/parquet", dest, revision=version, md5=md5s,
              matches_build2_pin=version == pins.FISHBASE_VERSION, database=name)
    return dest


USGS_COLUMNS = ["Species Name", "Scientific Name", "Minimum age at encounter", "Banding date", "Banding region",
                "Age at banding", "Sex", "Encounter date", "Encounter region", "How bird and band were obtained",
                "Condition of bird/band at encounter"]


def usgs_table(page: str) -> list[list[str]]:
    """The rows of the page's main table (id main_tbl), cell text with tags removed and entities decoded."""
    m = re.search(r'<table id="main_tbl".*?</table>', page, flags=re.S)
    if not m:
        raise Refused("the USGS longevity page has no table main_tbl")
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", m.group(0), flags=re.S):
        cells = [html.unescape(re.sub(r"<[^>]+>", "", c)).strip() for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", tr, flags=re.S)]
        cells = [re.sub(r"\s+", " ", c) for c in cells]
        if cells: rows.append(cells)
    return rows


def usgs(dest: Path, manifest: Path) -> Path:
    page_path = dest / "Longevity.php.html"
    fetch_to(pins.USGS_URL, page_path)
    rows = usgs_table(page_path.read_text(encoding="utf-8", errors="replace"))
    header, body = rows[0], rows[1:]
    if header != USGS_COLUMNS:
        raise Refused(f"the USGS longevity table has columns {header}, expected {USGS_COLUMNS}")
    out = dest / "longevity.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n"); w.writerow(header); w.writerows(body)
    h, c = md5_of(page_path), md5_of(out)
    _manifest(manifest, "usgs", pins.USGS_URL, dest, revision="unversioned page", md5={"page": h, "table_csv": c},
              records=len(body), matches_build2_pin={"page": h == pins.USGS_HTML_MD5, "table_csv": c == pins.USGS_CSV_MD5})
    if c != pins.USGS_CSV_MD5:
        print(f"WARNING: the USGS table differs from the one build 2 read (csv md5 {c} != {pins.USGS_CSV_MD5}); "
              f"{len(body)} records now")
    return dest


def euring(dest: Path, manifest: Path) -> Path:
    """The EURING longevity list of 2023-09-01, the Internet Archive's capture of the PDF (euring.org is in
    maintenance). Refused unless it has the md5 build 2 pinned."""
    dest.mkdir(parents=True, exist_ok=True)
    pdf = dest / pins.EURING_PDF
    fetch_to(pins.EURING_PDF_URL, pdf)
    h = md5_of(pdf)
    _manifest(manifest, "euring", pins.EURING_PDF_URL, dest, revision="list of 2023-09-01", md5={"pdf": h},
              matches_build2_pin=h == pins.EURING_MD5)
    if h != pins.EURING_MD5:
        raise Refused(f"the EURING PDF has md5 {h}, not build 2's {pins.EURING_MD5}")
    return pdf


# --- the old site's ids ------------------------------------------------------------------------
def legacy(dest: Path, manifest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    pages, page = 1, 1
    while page <= pages:
        p = dest / f"total_mtDNA_page{page:02d}.html"
        if not p.is_file():
            p.write_bytes(fetch(pins.LEGACY_URL.format(page=page)))
            time.sleep(0.5)
        pages = max([pages] + [int(x) for x in re.findall(r'href="\?page=(\d+)"', p.read_text(errors="replace"))])
        page += 1
    _manifest(manifest, "mitoage.info", pins.LEGACY_URL.format(page="N"), dest, revision="the old site (July 2015)", pages=pages)
    return dest


def taxdump_from_file(src: Path, dest_root: Path, manifest: Path) -> Path:
    """A kept new_taxdump.tar.gz of build 2's date, accepted only with its md5."""
    if not src.is_file():
        raise FileNotFoundError(f"no file at {src}")
    got = md5_of(src)
    if got != pins.TAXDUMP_BUILD2_MD5:
        raise Refused(f"{src} has md5 {got}; build 2's new_taxdump.tar.gz of {pins.TAXDUMP_BUILD2} is {pins.TAXDUMP_BUILD2_MD5}")
    dest = dest_root / f"current-{pins.TAXDUMP_BUILD2}"
    with tarfile.open(src) as tf:
        _extract(tf.getnames(), lambda n: tf.extractfile(n), dest)
    _manifest(manifest, "local-copy", str(src), dest, revision=f"new_taxdump of {pins.TAXDUMP_BUILD2}", md5=got,
              origin=pins.TAXDUMP_CURRENT_URL, matches_build2_pin=True)
    return dest


def contig(release: Path, folder: Path, manifest: Path) -> Path:
    """Complete metazoan NC_ records of a release that carry only a CONTIG line: their sequence by efetch."""
    from .refseq import completeness, is_complete, scan
    versions = [h.version for h, _ in scan(release) if not h.has_sequence and re.search(r"\bMetazoa\b", h.lineage)
                and is_complete(completeness(h))]
    folder.mkdir(parents=True, exist_ok=True)
    genbank(versions, folder, manifest, f"{len(versions)} CONTIG-only records of {release.parent.name}")
    return folder


def primary(species_file: Path, cited: Path, folder: Path, manifest: Path) -> Path:
    """The INSDC records named in the COMMENT ("derived from" / "identical to") of the cited records of the listed hagrids."""
    from .assemble import primary_accession
    from .build1 import load
    from .. import paths
    if not species_file.is_file():
        raise FileNotFoundError(f"no {species_file}: run `mitoage build2 primary-list` first")
    want = {int(x) for x in species_file.read_text().split()}
    raw = load(paths.B2_BUILD1)
    accs = []
    for s in raw.kept:
        if s.hagrid not in want: continue
        a = primary_accession((cited / f"{s.refseq}.gb").read_text())
        if a: accs.append(a)
    folder.mkdir(parents=True, exist_ok=True)
    genbank(accs, folder, manifest, f"INSDC source records of {len(accs)} cited RefSeq records")
    return folder
