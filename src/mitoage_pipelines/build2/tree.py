"""The dated species tree, as the site's pipeline builds it (site repository: pipeline/sources/tree.ts):

1. Names -> Open Tree of Life ids with TNRS (context Animals, exact matches only), in batches of 250; a
   species takes the id of its first name (current, then earlier names) not taken by another species.
2. The induced subtree of the synthetic tree for those ids (ids it does not know are dropped and the
   request repeated); a species mapped onto a broken taxon goes to the node OTL names for it.
3. Unary nodes collapsed, empty branches dropped, polytomies kept.
4. Each internal node dated by TimeTree 5 (MRCA of up to three NCBI taxids from each of its first 12 child
   clades; adjusted_age, or precomputed_age when that is 0).
5. Ages made monotone (a parent at least as old as its oldest child), undated nodes placed midway between
   their parent and their oldest child; branch length = parent age - node age.

Every answer is a JSON file in data/input/build2/tree, named as the site's cache names it (a hash of the
request), so a build reads the same answers as build 2 when it has them. `download tree-sources` asks
OTL and TimeTree for the answers that are missing; `build2 tree` never goes to the network.
"""
import hashlib
import json
import math
import sys
import time
import urllib.parse
from pathlib import Path

import polars as pl

from . import pins
from .download import Refused, fetch
from .lifespan import js_str


def h16(s: str) -> str:
    return hashlib.sha1(s.encode()).hexdigest()[:16]


class Answers:
    """The OTL and TimeTree answers in `folder`; `online` asks for missing ones (OTL 2/s, TimeTree 3/s)."""
    def __init__(self, folder: Path, online: bool):
        self.folder, self.online = folder, online
        self.fetched, self.failed, self.read = 0, 0, 0
        self.last = {"otl": 0.0, "tt": 0.0}
        self.used: list[str] = []          # every answer file the build read or fetched, in order

    def _wait(self, k, per_s):
        time.sleep(max(0.0, self.last[k] + 1 / per_s + 0.02 - time.monotonic())); self.last[k] = time.monotonic()

    def get(self, rel: str, fetcher, cache_failure=False):
        p = self.folder / rel
        self.used.append(rel)
        if p.is_file():
            self.read += 1
            return json.loads(p.read_text())
        if not self.online:
            raise FileNotFoundError(f"no {p}: the tree's answers are not all downloaded (`mitoage download tree-sources`)")
        v = fetcher()
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(v, separators=(",", ":")))
        self.fetched += 1
        return v

    def otl(self, path: str, body: dict):
        """POST to the OTL API; the JSON body is kept whatever the HTTP status (an unknown id is a 400 that names it)."""
        import urllib.error, urllib.request
        self._wait("otl", 2)
        req = urllib.request.Request(f"{pins.OTL_API}{path}", data=json.dumps(body).encode(),
                                     headers={"content-type": "application/json", "user-agent": "mitoage-pipelines"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r: raw = r.read()
        except urllib.error.HTTPError as e:
            raw = e.read()
        try:
            return json.loads(raw)
        except ValueError:
            raise Refused(f"OTL {path}: not JSON: {raw[:200]!r}")

    def timetree(self, ids: list[int]):
        self._wait("tt", 3)
        raw = fetch(f"{pins.TIMETREE_API}/mrca/id/{'+'.join(map(str, ids))}/json", tries=4)
        try:
            return json.loads(raw)
        except ValueError:
            return None          # a non-JSON answer is kept as "no age", as the site's code does


def parse_newick(s: str) -> dict:
    i = 0
    def label():
        nonlocal i
        out = ""
        if i < len(s) and s[i] == "'":
            i += 1
            while i < len(s) and s[i] != "'": out += s[i]; i += 1
            i += 1
        else:
            while i < len(s) and s[i] not in ",():;": out += s[i]; i += 1
        if i < len(s) and s[i] == ":":
            i += 1
            while i < len(s) and s[i] not in ",();": i += 1
        return out.strip()
    def node():
        nonlocal i
        kids = []
        if i < len(s) and s[i] == "(":
            i += 1; kids.append(node())
            while i < len(s) and s[i] == ",":
                i += 1; kids.append(node())
            i += 1
        return {"label": label(), "children": kids}
    return node()


def mrca_age(ans: Answers, taxids: list[int]):
    ids = sorted(set(taxids))
    if len(ids) < 2: return None
    try:
        res = ans.get(f"timetree/mrca_{h16('+'.join(map(str, ids)))}.json", lambda: ans.timetree(ids))
    except RuntimeError:
        ans.failed += 1
        return None              # network failures are not kept: the node stays undated
    if not res or len(res.get("found_ids") or []) < 2: return None
    for v in (res.get("adjusted_age"), res.get("precomputed_age")):
        try:
            age = float(v) if v is not None else None
        except (TypeError, ValueError):
            age = None
        if age is not None and math.isfinite(age) and age > 0: return age
    return None


def build(species: list[dict], ans: Answers) -> dict:
    """species: id, names (current first), taxid. Returns {version, root, stats}."""
    sys.setrecursionlimit(1_000_000)
    names = list(dict.fromkeys(n for s in species for n in s["names"]))
    ott: dict[str, int] = {}
    for i in range(0, len(names), 250):
        batch = names[i:i + 250]
        res = ans.get(f"tnrs_{h16('|'.join(batch))}.json", lambda: ans.otl("/tnrs/match_names",
                      {"names": batch, "context_name": "Animals", "do_approximate_matching": False}))
        for r in res["results"]:
            m = next((x for x in r["matches"] if not x["is_approximate_match"]), None)
            if m: ott[r["name"].lower()] = m["taxon"]["ott_id"]
    species_ott, ott_species = {}, {}
    for s in species:
        for n in s["names"]:
            o = ott.get(n.lower())
            if o is not None and o not in ott_species:
                species_ott[s["id"]] = o; ott_species[o] = s["id"]; break
    ids = sorted(species_ott.values())
    induced = None
    for _ in range(5):
        if induced: break
        res = ans.get(f"induced_{h16(','.join(map(str, ids)))}.json",
                      lambda: ans.otl("/tree_of_life/induced_subtree", {"ott_ids": ids, "label_format": "id"}))
        if "newick" in res: induced = res
        elif res.get("unknown"):
            bad = {int(k.replace("ott", "")) for k in res["unknown"]}
            ids = [x for x in ids if x not in bad]
        else: raise Refused("unexpected answer of OTL induced_subtree")
    if not induced: raise Refused("OTL induced subtree failed")
    about = ans.get("otl_about.json", lambda: ans.otl("/tree_of_life/about", {}))
    label_species = {f"ott{o}": sid for o, sid in ott_species.items()}
    for frm, to in (induced.get("broken") or {}).items():
        sid = label_species.get(frm)
        if sid is not None and to not in label_species: label_species[to] = sid
    taxid_of = {s["id"]: s["taxid"] for s in species}

    def convert(n):
        kids = [k for k in (convert(c) for c in n["children"]) if k is not None]
        sid = label_species.get(n["label"])
        if sid is not None:
            leaf = {"children": [], "leaf": sid, "age": 0}
            if not kids: return leaf
            kids.append(leaf)
        if not kids: return None
        if len(kids) == 1: return kids[0]
        return {"children": kids, "age": None}
    newick = induced["newick"].rstrip()
    root = convert(parse_newick(newick[:-1] if newick.endswith(";") else newick))
    if root is None: raise Refused("empty induced tree")

    def leaves(n):
        return [n["leaf"]] if "leaf" in n else [x for c in n["children"] for x in leaves(c)]
    internal = []
    def collect(n):
        for c in n["children"]: collect(c)
        if "leaf" not in n: internal.append(n)
    collect(root)
    dated = 0
    for k, n in enumerate(internal):
        reps = []
        for c in n["children"][:12]:
            reps += [t for t in (taxid_of.get(x) for x in leaves(c)) if isinstance(t, int)][:3]
        n["age"] = mrca_age(ans, reps)
        if n["age"] is not None: dated += 1
        if (k + 1) % 200 == 0: print(f"  TimeTree: {k + 1}/{len(internal)} nodes", flush=True)
    monotone, min_age = 0, {}
    def post(n):
        nonlocal monotone
        if "leaf" in n: return 0
        m = 0
        for c in n["children"]: m = max(m, post(c))
        min_age[id(n)] = m
        if n["age"] is not None and n["age"] < m:
            n["age"] = m; monotone += 1
        return n["age"] if n["age"] is not None else m
    post(root)
    interpolated = 0
    if root["age"] is None:
        root["age"] = min_age.get(id(root), 0) * 1.1 + 1; interpolated += 1
    def pre(n):
        nonlocal interpolated
        for c in n["children"]:
            if "leaf" not in c and c["age"] is None:
                c["age"] = ((n["age"] or 0) + min_age.get(id(c), 0)) / 2; interpolated += 1
            pre(c)
    pre(root)
    def phylo(n, parent_age):
        age = n["age"] or 0
        node = {"length": max(0, parent_age - age), "children": [phylo(c, age) for c in n["children"]]}
        if "leaf" in n: node["leaf"] = n["leaf"]
        return node
    tree = phylo(root, root["age"] or 0)
    n_leaves = len(leaves(root))
    return {"version": induced.get("synth_id") or about.get("synth_id") or "unknown", "root": tree,
            "stats": {"matchedNames": len(species_ott), "leaves": n_leaves, "internalNodes": len(internal), "datedByTimeTree": dated,
                      "interpolated": interpolated, "madeMonotone": monotone, "timetreeQueries": len(internal), "rootAge": root["age"] or 0}}


def to_json(v) -> str:
    """JSON as JavaScript's JSON.stringify writes it (numbers without a trailing .0)."""
    if isinstance(v, dict): return "{" + ",".join(f"{json.dumps(k)}:{to_json(x)}" for k, x in v.items()) + "}"
    if isinstance(v, list): return "[" + ",".join(to_json(x) for x in v) + "]"
    if isinstance(v, bool) or v is None: return json.dumps(v)
    if isinstance(v, (int, float)): return js_str(v)
    return json.dumps(v, ensure_ascii=False)


def species_input(species_parquet: Path) -> list[dict]:
    if not species_parquet.is_file():
        raise FileNotFoundError(f"no {species_parquet}: run `mitoage build2 assemble` first")
    sp = pl.read_parquet(species_parquet).sort("id")
    return [{"id": r["id"], "names": [r["scientific_name"]] + ([x for x in r["previous_scientific_names"].split("; ")] if r["previous_scientific_names"] else []),
             "taxid": r["ncbi_taxid"]} for r in sp.iter_rows(named=True)]


def main(species_parquet: Path, answers: Path, out: Path, online: bool) -> dict:
    sys.setrecursionlimit(1_000_000)
    ans = Answers(answers, online)
    t = build(species_input(species_parquet), ans)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(to_json({"version": t["version"], "root": t["root"]}))
    (out.parent / "phylogeny_stats.json").write_text(json.dumps({**t["stats"], "version": t["version"], "answers_read": ans.read,
                                                                 "answers_fetched": ans.fetched, "timetree_failures": ans.failed}, indent=1) + "\n")
    (out.parent / "phylogeny_answers.txt").write_text("".join(f"{u}\n" for u in dict.fromkeys(ans.used)))
    print(f"{t['version']}: {t['stats']}; answers read {ans.read}, fetched {ans.fetched}, TimeTree failures {ans.failed}")
    return t
