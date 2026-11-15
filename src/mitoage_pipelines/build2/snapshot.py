"""The snapshot dataset longevity-genie/mitoage-build2-sources on Hugging Face: copies of the inputs of
build 2 that no upstream archive keeps (RefSeq release 237, the taxonomy dump of 2026-11-13, the GenBank
records read by efetch, AnAge build 15, the USGS table, the OTL and TimeTree answers, the old site's
pages) and, for convenience, the FishBase/SeaLifeBase v26.06 tables read.

Reviewers only download (`mitoage download ... --from-snapshot`, no token needed). Building and uploading
it is a maintainer command (`mitoage snapshot build`, then `mitoage snapshot upload` with HF_TOKEN set in
the environment of that one command).
"""
import gzip
import hashlib
import json
import shutil
import tarfile
from datetime import datetime, timezone
from pathlib import Path

from . import pins
from .download import Refused, _manifest, md5_of

REPO = "longevity-genie/mitoage-build2-sources"
TAXDUMP_BUILD2_URL = pins.TAXDUMP_CURRENT_URL

LICENCES = {
    "ncbi": ("NCBI molecular data are free of use restrictions (NCBI website and data usage policies)", "NCBI RefSeq: O'Leary et al. 2016, Nucleic Acids Res 44:D733"),
    "taxonomy": ("NCBI Taxonomy is free of use restrictions (NCBI website and data usage policies)", "Schoch et al. 2020, Database baaa062"),
    "anage": (pins.ANAGE_LICENSE, "Tacutu et al. 2013, Nucleic Acids Res 41:D1027-D1033"),
    "fb": (pins.FISHBASE_LICENSE["fb"], "Froese R., Pauly D. (eds), FishBase, www.fishbase.org, version 06/2026"),
    "slb": (pins.FISHBASE_LICENSE["slb"], "Palomares M.L.D., Pauly D. (eds), SeaLifeBase, www.sealifebase.org, version 06/2026"),
    "usgs": (pins.USGS_LICENSE, "Bird Banding Laboratory, NABBP Longevity Records v2025.1, USGS EESC, Laurel MD"),
    "otl": ("CC0 1.0 (Open Tree of Life data)", "Hinchliff et al. 2015, PNAS 112:12764"),
    "timetree": ("Free public resource (TimeTree)", "Kumar et al. 2022, Mol Biol Evol 39:msac174"),
    "mitoage": (pins.BUILD1_LICENSE, "Toren et al. 2016, Nucleic Acids Res 44:D1262-D1265"),
}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()


def _gz_records(folder: Path, names: list[str] | None, dest: Path) -> int:
    """Records of a folder of <version>.gb files as one gzipped GenBank file (sorted by file name)."""
    files = sorted(folder.glob("*.gb")) if names is None else [folder / f"{n}.gb" for n in names]
    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(dest, "wb", compresslevel=9, mtime=0) as out:
        for f in files: out.write(f.read_bytes())
    return len(files)


def _tar(files: list[tuple[Path, str]], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as gz, tarfile.open(fileobj=gz, mode="w") as tf:
        for src, arc in files:
            info = tf.gettarinfo(str(src), arcname=arc); info.mtime = 0; info.uid = info.gid = 0; info.uname = info.gname = ""
            with open(src, "rb") as f: tf.addfile(info, f)


def build(staging: Path, b2_in: Path, sources: dict) -> list[dict]:
    """Copy the inputs into `staging` and describe each file. `sources`: where the copies come from
    (release237, taxdump_tar, tree_answers folder and the list of answers read, retrieval dates)."""
    if staging.exists(): shutil.rmtree(staging)
    staging.mkdir(parents=True)
    entries = []
    def add(rel, url, retrieved, lic, note, pin_md5=None):
        p = staging / rel
        md5 = md5_of(p)
        if pin_md5 and md5 != pin_md5: raise Refused(f"{rel}: md5 {md5} is not build 2's {pin_md5}")
        entries.append({"path": rel, "original_url": url, "retrieved_at": retrieved, "md5": md5, "sha256": sha256_of(p),
                        "bytes": p.stat().st_size, "licence": LICENCES[lic][0], "citation": LICENCES[lic][1], "note": note,
                        "matches_build2_pin": None if pin_md5 is None else True})
    rel = "refseq/release237/mitochondrion.1.genomic.gbff.gz"
    (staging / rel).parent.mkdir(parents=True); shutil.copyfile(sources["release237"], staging / rel)
    add(rel, pins.REFSEQ_URL, sources["release237_retrieved"], "ncbi",
        "NCBI RefSeq mitochondrion release 237 (file date 2026-09-03), the file build 2 read; NCBI keeps only the current release", pins.REFSEQ_MD5)
    rel = "taxonomy/new_taxdump_2026-11-13.tar.gz"
    (staging / rel).parent.mkdir(parents=True); shutil.copyfile(sources["taxdump_tar"], staging / rel)
    add(rel, TAXDUMP_BUILD2_URL, "2026-11-13T00:59:14Z", "taxonomy",
        "The current new_taxdump on 2026-11-13 (Last-Modified 2026-11-13T00:42:20Z), the dump build 2 read. NCBI's monthly archive keeps "
        "only the dumps of the first of each month; with the closest one (2026-11-01) the build has 2333 species instead of 2335",
        pins.TAXDUMP_BUILD2_MD5)
    for name, what, url_note in (("cited", "the 922 RefSeq records build 1 cites", "build 2 read the same text (byte-identical to its cached copies)"),
                                 ("contig", "complete records that carry only a CONTIG line in release 237 (6 of them are used)", "sequence by efetch rettype gbwithparts"),
                                 ("primary", "the INSDC records named in the COMMENT of the cited records the re-annotation test needs", "build2 primary-list chooses them")):
        rel = f"genbank/{name}.gbff.gz"
        n = _gz_records(b2_in / "genbank" / name, None, staging / rel)
        add(rel, pins.EFETCH_URL + "?db=nuccore&rettype=gbwithparts", sources["efetch_retrieved"], "ncbi", f"{n} records: {what}; {url_note}")
    rel = "anage/dataset.zip"
    (staging / rel).parent.mkdir(parents=True); shutil.copyfile(b2_in / "anage" / "dataset.zip", staging / rel)
    add(rel, pins.ANAGE_URL, sources["anage_retrieved"], "anage", "AnAge build 15 (2023-07-03) as downloaded; the URL is unversioned", pins.ANAGE_MD5)
    for db, folder in (("fb", "fishbase"), ("slb", "sealifebase")):
        for t in pins.FISHBASE_TABLES:
            rel = f"{folder}/{pins.FISHBASE_VERSION}/{t}.parquet"
            (staging / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(b2_in / folder / pins.FISHBASE_VERSION / f"{t}.parquet", staging / rel)
            add(rel, pins.SOURCECOOP_URL.format(db=db, version=pins.FISHBASE_VERSION, table=t), sources["fishbase_retrieved"], db,
                "mirror of a versioned source.coop file, for convenience; the downloader reads source.coop by default", pins.FISHBASE_MD5[db][t])
    for f, pin, note in (("longevity.csv", pins.USGS_CSV_MD5, "the table of the page as CSV (columns as on the page); identical to the table build 2 read"),
                         ("Longevity.php.html", None, "the page as retrieved (its HTML changes between retrievals; the table did not)")):
        rel = f"usgs_bbl/{f}"
        (staging / rel).parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(b2_in / "usgs_bbl" / f, staging / rel)
        add(rel, pins.USGS_URL, sources["usgs_retrieved"], "usgs", note, pin)
    rel = "tree/answers.tar.gz"
    _tar([(sources["tree_answers"] / u, u) for u in sources["tree_used"]], staging / rel)
    add(rel, f"{pins.OTL_API} (tnrs/match_names, tree_of_life/induced_subtree, tree_of_life/about); {pins.TIMETREE_API}/mrca/id/<ids>/json",
        sources["tree_retrieved"], "otl", f"the {len(sources['tree_used'])} answers the dated tree of build 2 was built from (OTL synthetic tree opentree16.1; "
        "TimeTree ages; file names are the hashes of the requests). Licences: OTL data CC0 1.0; TimeTree: free public resource, cite Kumar et al. 2022")
    rel = "legacy/total_mtDNA_pages.tar.gz"
    _tar([(p, p.name) for p in sorted((b2_in / "legacy").glob("*.html"))], staging / rel)
    add(rel, pins.LEGACY_URL.format(page="N"), sources["legacy_retrieved"], "mitoage", "the old MitoAge site's export pages, for its species ids (legacy_id)")
    (staging / "MANIFEST.json").write_text(json.dumps(entries, indent=1) + "\n")
    (staging / "README.md").write_text(card(entries))
    (staging / ".gitattributes").write_text("*.gz filter=lfs diff=lfs merge=lfs -text\n*.zip filter=lfs diff=lfs merge=lfs -text\n*.parquet filter=lfs diff=lfs merge=lfs -text\n")
    return entries


def card(entries: list[dict]) -> str:
    rows = "\n".join(f"| `{e['path']}` | {e['bytes']:,} | {e['retrieved_at']} | `{e['md5']}` | {e['original_url']} | {e['licence']} | {e['citation']} |"
                     for e in entries)
    notes = "\n".join(f"- `{e['path']}`: {e['note']}." for e in entries)
    return f"""---
license: other
license_name: mixed
license_link: LICENSE
pretty_name: MitoAge build 2 sources (snapshot)
tags:
- mitochondria
- longevity
- reproducibility
---

# MitoAge build 2 sources (snapshot)

A snapshot of the inputs of MitoAge build 2 (2026-11-13) that no upstream archive keeps, for reproducing it.
NCBI keeps only the current RefSeq release and monthly taxonomy dumps; AnAge, the USGS Bird Banding
Laboratory table and the TimeTree and Open Tree of Life answers are unversioned. Each file is described below
with its original URL, when it was retrieved, its checksums, its licence and how to cite it; `MANIFEST.json`
holds the same in machine-readable form.

- The build itself: [longevity-genie/mitoage](https://huggingface.co/datasets/longevity-genie/mitoage)
  (build 2 with the EURING bird records at revision `{pins.HF_BUILD2_REVISION}`, build 1 in its `build1/` folder).
- The pipeline that built it: [dna-seq/mitoage](https://github.com/dna-seq/mitoage) (TypeScript).
- The independent reproduction that reads this snapshot:
  [dna-seq/mitoage-pipelines](https://github.com/dna-seq/mitoage-pipelines),
  [docs/BUILD2_REPRODUCTION.md](https://github.com/dna-seq/mitoage-pipelines/blob/main/docs/BUILD2_REPRODUCTION.md).
  Its downloaders read the original sources by default and this snapshot with `--from-snapshot`; reading it
  needs no token.

Not here, because an upstream archive keeps it: the NCBI taxonomy monthly archive, the FishBase and
SeaLifeBase snapshots on source.coop (mirrored here only for convenience), the EURING longevity list of
2023-09-01 (the Internet Archive's capture of 2026-10-11 of `EURING_longevity_list_20230901.pdf`, md5
`{pins.EURING_MD5}`; euring.org shows a maintenance page, and the list states no licence), and the Hugging
Face build itself.

**Licences are per file.** FishBase and SeaLifeBase tables are CC BY-NC 4.0 (non-commercial use, with
attribution). AnAge is subject to the HAGR conditions of use (https://genomics.senescence.info/legal.html).
NCBI data have no use restrictions. The USGS table is a U.S. Government work. Open Tree of Life data are
CC0; TimeTree is a free public resource.

| file | bytes | retrieved (UTC) | md5 | original URL | licence | cite |
| --- | ---: | --- | --- | --- | --- | --- |
{rows}

## Notes per file

{notes}
"""


def upload(staging: Path, token_env: str = "HF_TOKEN") -> str:
    """Create the dataset (public) when needed and upload `staging` as one commit. The token comes from the
    environment variable HF_TOKEN of this command, or else from the machine's stored Hugging Face login, which
    huggingface_hub looks up itself; it is never read here, printed or written."""
    import os
    from huggingface_hub import HfApi
    api = HfApi(token=os.environ.get(token_env) or None)
    if api.whoami().get("name") is None:
        raise Refused(f"no Hugging Face login: set {token_env} or log in; uploading the snapshot is a maintainer step")
    api.create_repo(REPO, repo_type="dataset", private=False, exist_ok=True)
    info = api.upload_folder(repo_id=REPO, repo_type="dataset", folder_path=str(staging),
                             commit_message="MitoAge build 2 sources: snapshot of the inputs no upstream archive keeps")
    return info.oid


def fetch(path_in_repo: str, revision: str, dest: Path) -> tuple[Path, dict]:
    """One file of the snapshot (no token), checked against the snapshot's MANIFEST.json at the same revision."""
    from huggingface_hub import hf_hub_download
    man = json.loads(Path(hf_hub_download(REPO, "MANIFEST.json", repo_type="dataset", revision=revision, token=False)).read_text())
    entry = next((e for e in man if e["path"] == path_in_repo), None)
    if entry is None: raise Refused(f"{REPO}@{revision} has no {path_in_repo}")
    got = Path(hf_hub_download(REPO, path_in_repo, repo_type="dataset", revision=revision, token=False))
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(got, dest)
    md5 = md5_of(dest)
    if md5 != entry["md5"]: raise Refused(f"{REPO}@{revision}/{path_in_repo}: md5 {md5}, its manifest says {entry['md5']}")
    return dest, entry


def record(manifest: Path, folder: Path, path_in_repo: str, revision: str, entry: dict, **extra) -> dict:
    return _manifest(manifest, "huggingface-snapshot", f"datasets/{REPO}/{path_in_repo}", folder, revision=revision,
                     md5=entry["md5"], original_url=entry["original_url"], originally_retrieved=entry["retrieved_at"], **extra)


def unpack_records(gz: Path, folder: Path) -> int:
    """A gzipped multi-record GenBank file back into one <version>.gb file per record."""
    from .download import split_genbank, version_of
    folder.mkdir(parents=True, exist_ok=True)
    recs = split_genbank(gzip.open(gz, "rt").read())
    for r in recs: (folder / f"{version_of(r)}.gb").write_text(r)
    return len(recs)


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _untar(tar: Path, dest: Path) -> int:
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar) as tf:
        members = tf.getmembers()
        for m in members:
            if not (dest / m.name).resolve().is_relative_to(dest.resolve()) or not m.isfile():
                raise Refused(f"{tar}: unexpected member {m.name}")
        tf.extractall(dest, members=members, filter="data")
    return len(members)


def download(kind: str, revision: str, manifest: Path) -> None:
    """The snapshot's copy of one kind of input into the same place the upstream downloader writes it."""
    from .. import paths
    from .download import _extract
    tmp = paths.B2_IN / ".snapshot"
    if kind == "refseq":
        rel = "refseq/release237/mitochondrion.1.genomic.gbff.gz"
        dest, e = fetch(rel, revision, paths.refseq_release(pins.REFSEQ_RELEASE))
        if e["md5"] != pins.REFSEQ_MD5: raise Refused(f"{rel}: md5 {e['md5']} is not build 2's {pins.REFSEQ_MD5}")
        (dest.parent / "RELEASE_NUMBER").write_text(pins.REFSEQ_RELEASE + "\n")
        (dest.parent / "label.json").write_text(json.dumps({"label": pins.REFSEQ_VERSION}) + "\n")
        record(manifest, dest.parent, rel, revision, e, matches_build2_pin=True)
    elif kind == "taxonomy":
        rel = "taxonomy/new_taxdump_2026-11-13.tar.gz"
        tar, e = fetch(rel, revision, tmp / "new_taxdump.tar.gz")
        if e["md5"] != pins.TAXDUMP_BUILD2_MD5: raise Refused(f"{rel}: md5 {e['md5']} is not build 2's {pins.TAXDUMP_BUILD2_MD5}")
        dest = paths.taxdump(f"current-{pins.TAXDUMP_BUILD2}")
        with tarfile.open(tar) as tf: _extract(tf.getnames(), lambda n: tf.extractfile(n), dest)
        tar.unlink()
        record(manifest, dest, rel, revision, e, matches_build2_pin=True)
    elif kind == "lifespan":
        for rel, dest, pin in [("anage/dataset.zip", paths.B2_ANAGE, pins.ANAGE_MD5),
                               ("usgs_bbl/longevity.csv", paths.B2_USGS / "longevity.csv", pins.USGS_CSV_MD5),
                               ("usgs_bbl/Longevity.php.html", paths.B2_USGS / "Longevity.php.html", None)] + \
                              [(f"{f}/{pins.FISHBASE_VERSION}/{t}.parquet", (paths.B2_FISHBASE if db == "fb" else paths.B2_SEALIFEBASE) / pins.FISHBASE_VERSION / f"{t}.parquet",
                                pins.FISHBASE_MD5[db][t]) for db, f in (("fb", "fishbase"), ("slb", "sealifebase")) for t in pins.FISHBASE_TABLES]:
            d, e = fetch(rel, revision, dest)
            if pin and e["md5"] != pin: raise Refused(f"{rel}: md5 {e['md5']} is not build 2's {pin}")
            record(manifest, d.parent, rel, revision, e, matches_build2_pin=None if pin is None else True)
    elif kind in ("cited", "contig", "primary"):
        rel = f"genbank/{kind}.gbff.gz"
        gz, e = fetch(rel, revision, tmp / f"{kind}.gbff.gz")
        n = unpack_records(gz, paths.B2_GENBANK / kind); gz.unlink()
        record(manifest, paths.B2_GENBANK / kind, rel, revision, e, records=n)
    elif kind == "tree":
        rel = "tree/answers.tar.gz"
        tar, e = fetch(rel, revision, tmp / "answers.tar.gz")
        n = _untar(tar, paths.B2_TREE_IN); tar.unlink()
        record(manifest, paths.B2_TREE_IN, rel, revision, e, answers=n)
    elif kind == "legacy":
        rel = "legacy/total_mtDNA_pages.tar.gz"
        tar, e = fetch(rel, revision, tmp / "legacy.tar.gz")
        n = _untar(tar, paths.B2_LEGACY); tar.unlink()
        record(manifest, paths.B2_LEGACY, rel, revision, e, pages=n)
    else:
        raise ValueError(f"no snapshot of {kind}")
    if tmp.exists(): shutil.rmtree(tmp)
