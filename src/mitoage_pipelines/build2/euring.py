"""The EURING list of longevity records for European birds, list of 2023-09-01 (Fransson, Kolehmainen, Moss &
Robinson 2023), as build 2 reads it since the site repository's decision D15: the PDF the Internet Archive
captured on 2026-10-11 (euring.org shows a maintenance page), md5 pinned in pins.py.

The PDF is a spreadsheet exported by Excel: simple TrueType fonts in WinAnsi encoding, Flate-compressed content
streams, each table cell drawn at its own position. This module reads it with the standard library only and
refuses anything outside that (another filter, a composite font, a text operator it does not know).

The table: a bold row (common name, then the species' oldest record) and a row under it with the scientific
name (italic; one, Ciconia ciconia, is set upright) and the second-oldest record, when there is one. Columns:
country, ringing scheme, ring number, age ("37y 0m"), finding circumstances. Every age is a minimum.
"""
import re
import zlib
from dataclasses import dataclass
from pathlib import Path

from . import pins

# ---------------------------------------------------------------------------------------------- PDF
_OBJ = re.compile(rb"(\d+)\s+\d+\s+obj\b")
_REF = r"(\d+)\s+\d+\s+R"


@dataclass(frozen=True)
class Run:
    page: int
    x: float
    y: float
    font: str
    text: str


def _objects(data: bytes) -> dict[int, tuple[bytes, bytes | None]]:
    """Object number -> (dictionary text, raw stream bytes or None), objects in object streams included."""
    objs: dict[int, tuple[bytes, bytes | None]] = {}
    pos = 0
    while (m := _OBJ.search(data, pos)):
        end = data.find(b"endobj", m.end())
        if end < 0:
            raise ValueError(f"PDF object {m.group(1).decode()} has no endobj")
        body, stream = data[m.end():end], None
        s = re.search(rb"\bstream\r?\n", body)
        if s:
            stream = data[m.end() + s.end(): data.rfind(b"endstream", m.end(), end)]
            body = body[: s.start()]
        objs[int(m.group(1))] = (body, stream)
        pos = end
    for body, stream in list(objs.values()):
        if stream is None or not re.search(rb"/Type\s*/ObjStm\b", body):
            continue
        packed = _decode(body, stream, objs)
        n = int(re.search(rb"/N\s+(\d+)", body).group(1))
        first = int(re.search(rb"/First\s+(\d+)", body).group(1))
        head = [int(t) for t in packed[:first].split()]
        for k in range(n):
            start = first + head[2 * k + 1]
            stop = first + head[2 * k + 3] if k + 1 < n else len(packed)
            objs.setdefault(head[2 * k], (packed[start:stop], None))
    return objs


def _decode(body: bytes, stream: bytes, objs) -> bytes:
    ref = re.search(rb"/Length\s+" + _REF.encode(), body)
    length = int(objs[int(ref.group(1))][0].strip()) if ref else int(re.search(rb"/Length\s+(\d+)", body).group(1))
    raw = stream[:length]
    if re.search(rb"/Filter\s*\[", body):
        raise ValueError("a PDF stream with a filter array")
    f = re.search(rb"/Filter\s*/(\w+)", body)
    if f is None:
        return raw
    if f.group(1) == b"FlateDecode":
        return zlib.decompress(raw)
    raise ValueError(f"PDF stream filter {f.group(1).decode()} is not supported")


def _dict_of(body: bytes, key: str, objs) -> bytes | None:
    """The dictionary under /key, inline or by reference."""
    i = re.search(rb"/" + key.encode() + rb"\s*<<", body)
    if i:
        depth, j = 0, i.end() - 2
        while j < len(body):
            if body.startswith(b"<<", j): depth += 1; j += 2; continue
            if body.startswith(b">>", j):
                depth -= 1; j += 2
                if depth == 0: return body[i.end(): j - 2]
                continue
            j += 1
        raise ValueError(f"unbalanced /{key}")
    r = re.search(rb"/" + key.encode() + rb"\s+" + _REF.encode(), body)
    return objs[int(r.group(1))][0] if r else None


def _fonts(page: int, parents: dict[int, int], objs) -> dict[str, str]:
    """Font resource name -> BaseFont, with /Resources inherited from the page tree."""
    n: int | None = page
    while n is not None:
        res = _dict_of(objs[n][0], "Resources", objs)
        if res is not None:
            out = {}
            for name, ref in re.findall(rb"/(\w+)\s+" + _REF.encode(), _dict_of(res, "Font", objs) or b""):
                fb = objs[int(ref)][0]
                sub = re.search(rb"/Subtype\s*/(\w+)", fb)
                enc = re.search(rb"/Encoding\s*/(\w+)", fb)
                if sub is None or sub.group(1) not in (b"TrueType", b"Type1") or enc is None or enc.group(1) != b"WinAnsiEncoding":
                    raise ValueError(f"PDF font {name.decode()} is not a WinAnsi simple font")
                out[name.decode()] = re.search(rb"/BaseFont\s*/([^\s/<>\[\]]+)", fb).group(1).decode()
            return out
        n = parents.get(n)
    return {}


_TOKEN = re.compile(rb"\s*(?:(%[^\r\n]*)|(\()|(<<|>>)|(<[0-9A-Fa-f\s]*>)|(/[^\s/<>\[\]()%{}]*)|([\[\]])|"
                    rb"([+-]?(?:\d+\.?\d*|\.\d+))|([^\s/<>\[\]()%{}]+))")
_ESC = {b"n": 10, b"r": 13, b"t": 9, b"b": 8, b"f": 12, b"(": 40, b")": 41, b"\\": 92}


def _literal(s: bytes, i: int) -> tuple[bytes, int]:
    """A PDF literal string starting after its '(' at i: its bytes and the index after its ')'."""
    out, depth = bytearray(), 1
    while i < len(s):
        c = s[i:i + 1]
        if c == b"\\":
            e = s[i + 1:i + 2]
            if e in _ESC: out.append(_ESC[e]); i += 2
            elif e.isdigit() and e < b"8":
                o = re.match(rb"[0-7]{1,3}", s[i + 1:]).group(0); out.append(int(o, 8)); i += 1 + len(o)
            elif e in (b"\r", b"\n"): i += 3 if s[i + 1:i + 3] == b"\r\n" else 2
            else: out += e; i += 2
            continue
        if c == b"(": depth += 1
        elif c == b")":
            depth -= 1
            if depth == 0: return bytes(out), i + 1
        out += c; i += 1
    raise ValueError("unterminated PDF string")


def _runs(content: bytes, page: int, fonts: dict[str, str]) -> list[Run]:
    def text(b: bytes) -> str:
        return b.decode("cp1252")  # WinAnsiEncoding; bytes it leaves undefined raise
    runs, stack, font, x, y, leading, i = [], [], "", 0.0, 0.0, 0.0, 0
    while i < len(content):
        m = _TOKEN.match(content, i)
        if m is None or m.end() == i:
            if content[i:].strip() == b"": break
            raise ValueError(f"unreadable content at byte {i} of page {page}")
        i = m.end()
        com, lp, dd, hexs, name, br, num, op = m.groups()
        if com is not None or dd is not None: continue
        if lp is not None:
            b, i = _literal(content, i); stack.append(("str", b)); continue
        if hexs is not None:
            h = re.sub(rb"\s", b"", hexs[1:-1]); h += b"0" * (len(h) % 2); stack.append(("str", bytes.fromhex(h.decode()))); continue
        if name is not None: stack.append(("name", name[1:].decode())); continue
        if br is not None: stack.append((br.decode(), None)); continue
        if num is not None: stack.append(("num", float(num))); continue
        o = op.decode()
        nums = [v for k, v in stack if k == "num"]
        if o == "BT": x = y = 0.0
        elif o == "Tf":
            nm = [v for k, v in stack if k == "name"]
            if not nm or nm[-1] not in fonts: raise ValueError(f"unknown font on page {page}")
            font = fonts[nm[-1]]
        elif o == "Tm": x, y = nums[-2], nums[-1]
        elif o == "Td": x += nums[-2]; y += nums[-1]
        elif o == "TD": leading = -nums[-1]; x += nums[-2]; y += nums[-1]
        elif o == "TL": leading = nums[-1]
        elif o == "T*": y -= leading
        elif o == "Tj": runs.append(Run(page, x, y, font, text([v for k, v in stack if k == "str"][-1])))
        elif o == "TJ":
            open_at = max(j for j, (k, _) in enumerate(stack) if k == "[")
            runs.append(Run(page, x, y, font, text(b"".join(v for k, v in stack[open_at:] if k == "str"))))
        elif o in ("'", '"'): raise ValueError(f"text operator {o} is not supported")
        stack = []
    return runs


def text_runs(data: bytes) -> list[Run]:
    """Every text run, page by page in page-tree order."""
    objs = _objects(data)
    root = int(re.search(rb"/Root\s+" + _REF.encode(), data[data.rfind(b"trailer"):]).group(1))
    pages, parents = [], {}

    def walk(n: int) -> None:
        body = objs[n][0]
        if re.search(rb"/Type\s*/Pages\b", body):
            for k in re.findall(_REF.encode(), re.search(rb"/Kids\s*\[([^\]]*)\]", body).group(1)):
                parents[int(k)] = n; walk(int(k))
        else:
            pages.append(n)
    walk(int(re.search(rb"/Pages\s+" + _REF.encode(), objs[root][0]).group(1)))
    out = []
    for p, n in enumerate(pages, start=1):
        body = objs[n][0]
        arr = re.search(rb"/Contents\s*\[([^\]]*)\]", body)
        refs = re.findall(_REF.encode(), arr.group(1)) if arr else [re.search(rb"/Contents\s+" + _REF.encode(), body).group(1)]
        content = b"\n".join(_decode(*objs[int(r)], objs) for r in refs)
        out += _runs(content, p, _fonts(n, parents, objs))
    return out


# ---------------------------------------------------------------------------------------------- the list
@dataclass(frozen=True)
class Record:
    common_name: str
    scientific_name: str
    rank: int              # 1: the species' oldest record, 2: the second oldest
    country: str
    scheme: str
    ring_no: str
    age_text: str
    years: float | None
    age_note: str | None
    finding: str


def collapse(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def age(s: str) -> tuple[float | None, str | None]:
    """'37y 0m' -> 37.0; a month that is not a number below 12 keeps the whole years, and says so."""
    m = re.match(r"^\s*(\d+)\s*y\s*(\S+)\s*$", s)
    if not m:
        return None, f'age "{s.strip()}" not readable'
    mm = re.fullmatch(r"(\d+)m", m.group(2))
    if mm and int(mm.group(1)) < 12:
        return int(m.group(1)) + int(mm.group(1)) / 12, None
    return float(int(m.group(1))), f'months "{m.group(2)}" not readable; whole years used'


def _column(x: float) -> str:
    # cells are centred, so their x varies: bands between the columns of the table
    return "country" if x < 240 else "scheme" if x < 290 else "ring_no" if x < 340 else "age_text" if x < 390 else "finding"


def parse(runs: list[Run]) -> list[Record]:
    """The records of the list, top to bottom. Refuses a row it cannot place."""
    rows: list[list[Run]] = []
    for r in sorted(runs, key=lambda r: (r.page, -r.y, r.x)):
        if rows and rows[-1][0].page == r.page and abs(rows[-1][0].y - r.y) < 1: rows[-1].append(r)
        else: rows.append([r])
    is_header = lambda row: row[0].text.startswith("Species / Scientific name")
    start = next((i for i, row in enumerate(rows) if is_header(row)), None)
    if start is None:
        raise ValueError("EURING list: no table header")

    def cells(row: list[Run]) -> dict:
        c = dict.fromkeys(("country", "scheme", "ring_no", "age_text", "finding"), "")
        for r in sorted(row, key=lambda r: r.x)[1:]:
            k = _column(r.x)
            if c[k]: raise ValueError(f"EURING list: two {k} cells on page {r.page}")
            c[k] = collapse(r.text)
        if not c["age_text"]: raise ValueError(f"EURING list: '{row[0].text}' on page {row[0].page} has no age")
        c["years"], c["age_note"] = age(c["age_text"])
        return c

    out: list[Record] = []
    common = first_cells = None
    for row in rows[start + 1:]:
        if is_header(row): continue
        row = sorted(row, key=lambda r: r.x)
        lead = row[0]
        if lead.x > 30: raise ValueError(f"EURING list: a row without a name on page {lead.page}")
        name = collapse(lead.text)
        if "Bold" in lead.font:
            if common is not None: raise ValueError(f"EURING list: {common} has no scientific name")
            common, first_cells = name, cells(row)
        else:
            if common is None: raise ValueError(f"EURING list: '{name}' follows no common name")
            out.append(Record(common, name, 1, **first_cells))
            if len(row) > 1: out.append(Record(common, name, 2, **cells(row)))
            common = None
    if common is not None: raise ValueError(f"EURING list: {common} has no scientific name")
    return out


def load(pdf: Path):
    """The list as observations, after checking the file is the pinned one."""
    from .download import md5_of
    from .lifespan import Obs
    if not pdf.is_file():
        raise FileNotFoundError(f"no {pdf}: run `mitoage download lifespan` first")
    if md5_of(pdf) != pins.EURING_MD5:
        raise ValueError(f"{pdf} has md5 {md5_of(pdf)}, not the pinned {pins.EURING_MD5}")
    obs = []
    for r in parse(text_runs(pdf.read_bytes())):
        note = f" ({r.age_note})" if r.age_note else ""
        obs.append(Obs("EURING", pins.EURING_VERSION, f"euring20230901:{r.scheme}:{r.ring_no}", r.scientific_name,
                       r.common_name, r.years, "wild", "acceptable",
                       f"single individual ({'oldest' if r.rank == 1 else 'second oldest'} ringed in Europe; minimum age)", None,
                       f"age {r.age_text}{note}; {r.country}, scheme {r.scheme}, ring {r.ring_no}; {r.finding}",
                       pins.EURING_REFERENCE, pins.EURING_LICENSE))
    return obs
