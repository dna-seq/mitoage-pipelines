"""Hibernation of the build's species, from two pinned sources (H23, H24, H25).

- Bats: the trait table of Wilkinson & Adams 2019 ("Recurrent evolution of extreme longevity in bats", Biol Lett
  15:20180860, doi:10.1098/rsbl.2018.0860), deposited in the University of Maryland's repository DRUM
  (hdl:1903/21501) as one tab-separated text file. Only its `Genus_species` and `Hibernation` columns are read:
  the file also carries longevity, which the tests of H24 and H25 take from the build, never from here.
- Mammals: the supplementary table of Heldstab et al. 2018 ("Hibernation constrains brain size evolution in
  mammals", J Evol Biol, doi:10.1111/jeb.13353), Dryad doi:10.5061/dryad.753d06g, CC0. Dryad serves files only to
  a browser or a token holder, so the file is downloaded by hand from the dataset page and placed at
  `MAMMALS_FILE`; `fetch` accepts it only with the md5 Dryad's own file listing gives. One sheet, "Appendix S1",
  1,104 species; only `Genus species` and `hibernation` (1 or 0) are read. It carries no lifespan.

Matching to the build: a table name equal to the build's scientific name; otherwise a name that NCBI's taxonomy
(names.dmp, any name class) gives to the build species' own taxid. A build species reached by two table rows with
different states, or a table name reaching two build species, is "ambiguous" and its state unknown. A build species
not in the table is unknown, never "does not hibernate".
"""
import hashlib
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

BATS_URL = "https://drum.lib.umd.edu/bitstreams/fea60378-73ce-43ad-918f-9674cbff4dbe/download"
BATS_HANDLE = "http://hdl.handle.net/1903/21501"
BATS_FILE = "wilkinson_adams_2019_bats.txt"
BATS_SHA256 = "c2b402068f93991c58eb81527946c21601db15f482793c7913d050029441c2c0"   # as served on 2026-11-16

MAMMALS_PAGE = "https://datadryad.org/dataset/doi:10.5061/dryad.753d06g"
MAMMALS_FILE = "heldstab_2018_supplementary_material.xlsx"     # "Supplementary material.xlsx" on Dryad, renamed
MAMMALS_MD5 = "c5e58d7f194bc830bc7ba854ad54aabe"                # Dryad's file listing (API v2, version 9801)

STATES = {"Yes": 1, "No": 0}


class Refused(RuntimeError):
    """A source that is not what it should be."""


def digest(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_bats(path: Path) -> pd.DataFrame:
    """name (genus and species with a space), family, hibernation (1/0). The file has CR line ends and tabs.
    Refuses a value of Hibernation other than Yes or No, an empty name, and a name given twice."""
    text = path.read_bytes().decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    rows = [line.split("\t") for line in text.split("\n") if line.strip()]
    head, body = rows[0], rows[1:]
    try:
        i_name, i_fam, i_hib = head.index("Genus_species"), head.index("Family"), head.index("Hibernation")
    except ValueError as e:
        raise Refused(f"{path.name}: no column {e}") from e
    out = []
    for n, r in enumerate(body, start=2):
        name, value = r[i_name].strip().replace("_", " "), r[i_hib].strip()
        if not name:
            raise Refused(f"{path.name} line {n}: empty Genus_species")
        if value not in STATES:
            raise Refused(f"{path.name} line {n}: Hibernation {value!r} is neither Yes nor No")
        out.append((name, r[i_fam].strip(), STATES[value]))
    t = pd.DataFrame(out, columns=["name", "family", "hibernation"])
    dup = t.name[t.name.duplicated()].tolist()
    if dup:
        raise Refused(f"{path.name}: names given twice: {', '.join(sorted(set(dup)))}")
    return t


XLSX_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def xlsx_rows(path: Path, sheet: str = "xl/worksheets/sheet1.xml") -> list[list[str]]:
    """The cells of one sheet as text, rows in order, each row as wide as its last filled cell; an empty cell is "".
    Shared and inline strings are resolved; numbers are kept as the file writes them. (No spreadsheet package is a
    dependency of this repository, and this file needs nothing more.)"""
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", XLSX_NS):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{XLSX_NS['m']}}}t")))
        root = ET.fromstring(z.read(sheet))
    out = []
    for row in root.iter(f"{{{XLSX_NS['m']}}}row"):
        cells = {}
        for c in row.findall("m:c", XLSX_NS):
            col = 0
            for ch in re.match(r"[A-Z]+", c.get("r")).group():
                col = col * 26 + ord(ch) - 64
            kind, v = c.get("t"), c.find("m:v", XLSX_NS)
            if kind == "s":
                text = shared[int(v.text)]
            elif kind == "inlineStr":
                text = "".join(t.text or "" for t in c.iter(f"{{{XLSX_NS['m']}}}t"))
            else:
                text = v.text if v is not None and v.text is not None else ""
            cells[col - 1] = text
        out.append([cells.get(i, "") for i in range(max(cells) + 1)] if cells else [])
    return out


def read_mammals(path: Path) -> pd.DataFrame:
    """name, group (the order), hibernation (1/0) from Heldstab et al.'s sheet. Refuses a hibernation value other than
    1 or 0, an empty name, and a name given twice."""
    rows = xlsx_rows(path)
    head = rows[0]
    try:
        i_name, i_grp, i_hib = head.index("Genus species"), head.index("Group"), head.index("hibernation")
    except ValueError as e:
        raise Refused(f"{path.name}: no column {e}") from e
    out = []
    for n, r in enumerate(rows[1:], start=2):
        r = r + [""] * (len(head) - len(r))
        name, value = " ".join(r[i_name].split()), r[i_hib].strip()
        if not name:
            raise Refused(f"{path.name} row {n}: empty Genus species")
        if value not in ("0", "1"):
            raise Refused(f"{path.name} row {n}: hibernation {value!r} is neither 1 nor 0")
        out.append((name, r[i_grp].strip(), int(value)))
    t = pd.DataFrame(out, columns=["name", "group", "hibernation"])
    dup = t.name[t.name.duplicated()].tolist()
    if dup:
        raise Refused(f"{path.name}: names given twice: {', '.join(sorted(set(dup)))}")
    return t


def combine(species: pd.DataFrame, mammals: pd.DataFrame | None, bats: pd.DataFrame | None,
            names_dmp: Path | None) -> pd.DataFrame:
    """One state per build species from both sources (each matched by `match`): slug, state, how, and each source's
    own state. Known in one source only: that state. Known in both and equal: that state. Known in both and different:
    unknown, "sources disagree". Known in neither: unknown."""
    cols = {}
    for key, table in (("mammals", mammals), ("bats", bats)):
        if table is not None:
            cols[key] = match(species, table, names_dmp).set_index("slug").state
    out = pd.DataFrame(cols, index=pd.Index(species.slug, name="slug")).astype("Int64")
    known = out.notna().sum(axis=1)
    state, how = pd.Series(pd.NA, index=out.index, dtype="Int64"), pd.Series("in neither", index=out.index)
    one = known == 1
    state[one] = out[one].bfill(axis=1).iloc[:, 0]
    how[one] = "one source"
    two = known == 2
    agree = two & out.nunique(axis=1).eq(1)
    state[agree] = out[agree].iloc[:, 0]
    how[agree] = "both agree"
    how[two & ~agree] = "sources disagree"
    return out.assign(state=state, how=how).reset_index()


def ncbi_taxids(names_dmp: Path, wanted: set[str]) -> dict[str, set[int]]:
    """name -> the taxids NCBI gives that name under any name class, for the names in `wanted` only."""
    out: dict[str, set[int]] = {}
    with open(names_dmp, encoding="utf-8") as f:
        for line in f:
            parts = line.split("\t|\t")
            if len(parts) > 1 and parts[1] in wanted:
                out.setdefault(parts[1], set()).add(int(parts[0]))
    return out


def match(species: pd.DataFrame, table: pd.DataFrame, names_dmp: Path | None) -> pd.DataFrame:
    """One row per build species (columns slug, scientific_name, ncbi_taxid): slug, state (Int64, NA = unknown),
    how ("exact name", "NCBI synonym", "not in the table", "ambiguous ...") and the table name used."""
    by_name = dict(zip(table.name, table.hibernation))
    taxids = ncbi_taxids(names_dmp, set(table.name)) if names_dmp is not None else {}
    by_taxid: dict[int, list[str]] = {}
    for name, tids in taxids.items():
        for t in tids:
            by_taxid.setdefault(t, []).append(name)
    exact = set(species.scientific_name)
    rows = []
    for slug, sci, tid in zip(species.slug, species.scientific_name, species.ncbi_taxid):
        if sci in by_name:
            rows.append((slug, by_name[sci], "exact name", sci))
            continue
        cands = [] if pd.isna(tid) else [n for n in by_taxid.get(int(tid), []) if n not in exact]
        states = {by_name[n] for n in cands}
        if not cands:
            rows.append((slug, pd.NA, "not in the table", ""))
        elif len(states) > 1:
            rows.append((slug, pd.NA, "ambiguous: table rows disagree", "; ".join(sorted(cands))))
        else:
            rows.append((slug, states.pop(), "NCBI synonym", "; ".join(sorted(cands))))
    out = pd.DataFrame(rows, columns=["slug", "state", "how", "table_name"]).astype({"state": "Int64"})
    used = out.table_name[out.how == "NCBI synonym"]
    twice = used[used.duplicated(keep=False)]
    out.loc[out.table_name.isin(twice), ["state", "how"]] = [pd.NA, "ambiguous: one table name, two build species"]
    return out


def fetch(dest: Path, manifest: Path) -> dict[str, Path | None]:
    """The bat table from DRUM, refused unless its SHA-256 is the pinned one; and the mammal table if it has been
    placed by hand, refused unless its md5 is Dryad's. A missing mammal table is said, not hidden."""
    from ..build2 import download as dl
    dest.mkdir(parents=True, exist_ok=True)
    bats = dest / BATS_FILE
    dl.fetch_to(BATS_URL, bats)
    got = digest(bats, "sha256")
    if got != BATS_SHA256:
        raise Refused(f"{BATS_URL}: SHA-256 {got}, pinned {BATS_SHA256}")
    dl._manifest(manifest, "drum", BATS_HANDLE, dest, revision="as served 2026-11-16", sha256={BATS_FILE: got},
                 note="Wilkinson & Adams 2019 trait table; only Genus_species and Hibernation are read")
    mammals = dest / MAMMALS_FILE
    if not mammals.is_file():
        print(f"not placed: {mammals}. Download 'Supplementary material.xlsx' from {MAMMALS_PAGE} and save it there; "
              "until then the mammal-wide hibernation readouts are not assessed.")
        return {"bats": bats, "mammals": None}
    got = digest(mammals, "md5")
    if got != MAMMALS_MD5:
        raise Refused(f"{mammals}: md5 {got}, Dryad lists {MAMMALS_MD5}")
    dl._manifest(manifest, "dryad (placed by hand)", "doi:10.5061/dryad.753d06g", dest, revision="version 9801",
                 md5={MAMMALS_FILE: got})
    return {"bats": bats, "mammals": mammals}
