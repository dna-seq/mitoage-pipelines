"""Build 2 from build 1 and the sources: the status of the 922 build 1 species, the new species, the genome
rows recomputed with the build-2 rules, the changelog and the row provenance (site repository:
pipeline/build.ts, pipeline/steps/build2.ts; decisions in plans/DATA_BUILD2_DECISIONS.md and
docs/BUILD2_OPEN_QUESTIONS.md). Each stage is a function; `mitoage build2 <stage>` runs the stages up to it.
"""
import json
import re
import time
import unicodedata
from dataclasses import dataclass, field
from pathlib import Path

from . import build1 as b1mod
from . import diff as dmod
from . import genome as gm
from . import lifespan as ls
from . import pins
from .refseq import classify, one_per_species, scan, successor_for

NEW_ID_OFFSET = 10_000_000
POSSIBLE_DUPLICATES = [   # from build 2's audit, kept and flagged (open question h); not detected automatically
    ("Gadus ogac", "Gadus macrocephalus", "FishBase treats these NCBI species as one (Gadus macrocephalus); each has its own RefSeq mitogenome"),
    ("Panopea generosa", "Panopea abrupta", "both 168 y (AnAge vs SeaLifeBase); P. abrupta is widely treated as a synonym of P. generosa"),
    ("Euthynnus affinis", "Euthynnus alletteratus", "same genus, same maximum longevity (10 y) and the same AnAge common name"),
]
RENAME_CLASSES = {"genus_transfer", "epithet_spelling_or_gender", "genus_transfer_with_epithet_agreement"}
SRC_ORDER = ["AnAge", "FishBase", "SeaLifeBase", "EURING", "USGS_BBL", "MitoAge_build1"]


def mint_id(taxid: int) -> int:
    return NEW_ID_OFFSET + taxid


def slugify(s: str) -> str:
    s = "".join(c for c in unicodedata.normalize("NFD", s) if not 0x300 <= ord(c) <= 0x36F).lower()
    return re.sub(r"^-+|-+$", "", re.sub(r"[^a-z0-9]+", "-", s))


def name_change(b1: str, now: str | None) -> str | None:
    """How a later name relates to a build 1 name."""
    if not now or b1 == now: return None
    p1, p2 = b1.split(" "), now.split(" ")
    g1, e1 = p1[0], (p1[1] if len(p1) > 1 else None)
    g2, e2 = p2[0], (p2[1] if len(p2) > 1 else None)
    if g1 == g2 and e1 == e2: return None
    sim = (e1 or "")[:4] == (e2 or "")[:4]
    if g1 != g2 and e1 == e2: return "genus_transfer"
    if g1 == g2 and sim: return "epithet_spelling_or_gender"
    if g1 != g2 and sim: return "genus_transfer_with_epithet_agreement"
    return "different_species_name"


def primary_accession(text: str) -> str | None:
    """The INSDC record a RefSeq record's COMMENT names ("derived from X" / "identical to X")."""
    m = re.search(r"\nCOMMENT(.*?)\nFEATURES", text, flags=re.S)
    c = re.sub(r"\s+", " ", m.group(1)) if m else ""
    m = re.search(r"(?:identical to|derived from)\s+([A-Z]{1,2}_?\d+(?:\.\d+)?)", c)
    return m.group(1) if m else None


@dataclass
class Inputs:
    build1: Path
    release: Path
    release_label: str
    taxdump: Path
    taxonomy_label: str
    anage: Path
    fishbase: Path
    sealifebase: Path
    usgs: Path
    euring: Path
    genbank: Path               # cited/, contig/, primary/
    legacy: Path
    build_date: str = pins.BUILD_DATE


@dataclass
class State:
    """What the stages compute, in order."""
    inp: Inputs
    raw: b1mod.Build1 = None
    legacy: dict = field(default_factory=dict)
    cited: dict = field(default_factory=dict)          # hagrid -> Genome of the record build 1 cites
    cited_text: dict = field(default_factory=dict)
    anage: list = None
    anage_build: str = ""
    anage_match: dict = field(default_factory=dict)    # hagrid -> (entry or None, how)
    base: dict = field(default_factory=dict)           # hagrid -> the build 1 view of the species
    headers: list = None
    tax: object = None
    obs: list = None
    synonyms: dict = None
    matches: dict = None
    records: list = None
    log: list = field(default_factory=list)

    def step(self, msg: str, t0=[time.time()]):
        print(f"[{time.time() - t0[0]:6.1f}s] {msg}", flush=True)


# --- stage 1: build 1, its cited records, AnAge links --------------------------------------------
def load_build1(st: State) -> None:
    inp = st.inp
    st.step("build 1 (Hugging Face build1/) and the records it cites")
    st.raw = b1mod.load(inp.build1)
    for p in sorted(inp.legacy.glob("total_mtDNA_page*.html")):
        for m in re.finditer(r'<a href="/view/(\d+)/species/">([^<]+)</a>', p.read_text(errors="replace")):
            st.legacy[m.group(2).strip()] = int(m.group(1))
    if not st.legacy:
        print("  the old site's pages are not downloaded (`mitoage download legacy`): legacy_id stays empty")
    missing = []
    for s in st.raw.kept:
        p = inp.genbank / "cited" / f"{s.refseq}.gb"
        if not p.is_file(): missing.append(s.refseq); continue
        text = p.read_text()
        g = gm.annotate(text, gm.fallback_table(s.phylum), s.phylum == "Chordata")
        if g.record.accession_version != s.refseq: raise ValueError(f"{p} holds {g.record.accession_version}")
        st.cited[s.hagrid] = g; st.cited_text[s.hagrid] = text
    if missing:
        raise FileNotFoundError(f"{len(missing)} records cited by build 1 are not downloaded (`mitoage download cited`): {', '.join(missing[:5])}")
    st.anage = ls.load_anage(inp.anage)
    import zipfile
    with zipfile.ZipFile(inp.anage) as z:
        rel = next((n for n in z.namelist() if n.endswith("release.html")), None)
        txt = re.sub(r"<[^>]*>", " ", z.read(rel).decode("utf-8", "replace")) if rel else ""
    m = re.search(r"AnAge\s+(Build\s+\d+)", txt, flags=re.I)
    st.anage_build = m.group(1) if m else "unknown build"
    by_name = {e.scientific_name.lower(): e for e in st.anage}
    by_common: dict[str, list] = {}
    for e in st.anage: by_common.setdefault(e.common_name.lower(), []).append(e)
    for s in st.raw.kept:
        e = by_name.get(s.species.lower())
        if e: st.anage_match[s.hagrid] = (e, "scientific name")
        else:
            c = by_common.get(s.common_name.lower(), [])
            st.anage_match[s.hagrid] = (c[0], "common name") if len(c) == 1 else (None, "none")
    build1_view(st)


def build1_view(st: State) -> None:
    """Build 1 as build 2 starts from it: CSV values with the base-label fix, the cited record's features."""
    raw = st.raw
    used, slugs = set(), {}
    for s in sorted(raw.kept, key=lambda x: x.hagrid):
        slug = slugify(s.species)
        if slug in used: slug = f"{slug}-{slugify(s.refseq)}"
        used.add(slug); slugs[s.hagrid] = slug
    total_rows, pc_rows = raw.composition["total_mtDNA"], raw.composition.get("total_pc", {})
    for s in raw.kept:
        total = total_rows[s.hagrid][0]
        fx = b1mod.fixed_counts(total)
        e, _ = st.anage_match[s.hagrid]
        g = st.cited[s.hagrid]
        traits = ({**e.traits, "source": f"AnAge {st.anage_build} (HAGRID {e.hagrid})"} if e else
                  {**{k: None for k in ls.TRAITS}, "source": None})
        st.base[s.hagrid] = {
            "slug": slugs[s.hagrid], "gc_percent": gm.js_round(100 * (fx["g"] + fx["c"]) / (fx["a"] + fx["c"] + fx["g"] + fx["t"]), 3),
            "non_coding": total.non_coding if s.hagrid in pc_rows else gm.non_coding(g), "mtdna_size": total.length,
            "traits": traits, "legacy_id": st.legacy.get(s.species),
            "features": [{"name": f.name, "type": f.type, "start": f.start, "end": f.end, "strand": f.strand} for f in g.features],
        }


def build1_rows(st: State, hagrid: int) -> tuple[list, list]:
    """Build 1's own region and codon rows of one species (several CSV rows of a region: the one the cited
    record reproduces, else the closest in length)."""
    g = st.cited[hagrid]
    stats, usage = [], []
    for region in gm.REGIONS:
        rows = st.raw.composition.get(region, {}).get(hagrid, [])
        if not rows: continue
        chosen = rows[0]
        if len(rows) > 1:
            comp = gm.compute_region(g, region)
            if comp:
                match = next((r for r in rows if dmod.same_acgt(b1mod.fixed_counts(r), comp[0])), None)
                closest = sorted(rows, key=lambda r: abs(r.length - comp[1]))[0]
                chosen = match or closest
        stats.append(gm.make_stats(hagrid, region, b1mod.fixed_counts(chosen), "mitoage"))
    for region in gm.CODON_REGIONS:
        row = st.raw.codons.get(region, {}).get(hagrid)
        if row is None: continue
        counts = {c: row.get(c, 0) for c in gm.CODONS}
        usage.append(gm.usage_row(hagrid, region, counts, g.transl_table, "mitoage"))
    return stats, usage


# --- stage 2: the release, taxonomy and lifespan sources, names matched -------------------------
def load_sources(st: State) -> None:
    inp = st.inp
    st.step(f"RefSeq release headers ({inp.release_label})")
    st.headers = [h for h, _ in scan(inp.release)]
    print(f"  {len(st.headers)} records")
    st.step(f"NCBI Taxonomy ({inp.taxonomy_label})")
    from .taxonomy import Taxonomy
    tax = st.tax = Taxonomy(inp.taxdump)
    lineage_for = set()
    for h in st.headers:
        if h.taxid is None: continue
        sp = tax.species_of(h.taxid)
        if sp is not None: lineage_for.add(sp)
    for g in st.cited.values():
        if g.taxid is not None: lineage_for.add(tax.species_of(g.taxid) or tax.current(g.taxid))
    st.step("lifespan sources")
    st.anage, st.obs, st.synonyms = st.anage, *ls.load_all({"anage": inp.anage, "fishbase": inp.fishbase,
                                                             "sealifebase": inp.sealifebase, "usgs": inp.usgs,
                                                             "euring": inp.euring / pins.EURING_PDF})[1:]
    counts = {}
    for o in st.obs: counts[o.source] = counts.get(o.source, 0) + 1
    print(f"  observations: {counts}")
    name_taxids = set(lineage_for)
    for g in st.cited.values():
        if g.taxid is not None: name_taxids.add(tax.current(g.taxid))
    st.step("matching source names to NCBI species")
    tax.load_names(ls.wanted_name_keys(st.obs, st.synonyms), name_taxids)
    st.matches = ls.match_names(st.obs, tax, st.synonyms)
    by = {}
    for m in st.matches.values(): by[f"{m.source} {m.status}"] = by.get(f"{m.source} {m.status}", 0) + 1
    print("  " + ", ".join(f"{k}: {v}" for k, v in sorted(by.items())))
    st.records = classify(st.headers, tax)


# --- stage 3: status of the build 1 species and the candidates ----------------------------------
def plan(st: State) -> dict:
    tax, review = st.tax, []
    RS = f"NCBI RefSeq {st.inp.release_label}"
    by_taxid: dict[int, list] = {}
    for o in st.obs:
        m = st.matches.get(ls.rec_key(o.source, o.source_name, o.speccode))
        if ls.usable(m): by_taxid.setdefault(m.species_taxid, []).append(ls.matched(o, m))
    choice = one_per_species(st.records)
    record_by_acc = {r.header.accession: r for r in st.records}
    by_version = {h.version: h for h in st.headers}
    by_accession = {h.accession: h for h in st.headers}

    origs = []
    for s in st.raw.kept:
        g1 = st.cited[s.hagrid]
        rec_taxid = None if g1.taxid is None else tax.current(g1.taxid)
        sp_taxid = None if rec_taxid is None else tax.species_of(rec_taxid)
        entry = st.anage_match[s.hagrid][0]
        am = st.matches.get(ls.rec_key("AnAge", entry.scientific_name, None)) if entry else None
        origs.append({"s": s, "id": s.hagrid, "g1": g1, "rec_taxid": rec_taxid, "sp_taxid": sp_taxid, "eff": sp_taxid,
                      "entry": entry, "anage_taxid": am.species_taxid if am else None})
    sp_count = {}
    for o in origs:
        if o["sp_taxid"] is not None: sp_count[o["sp_taxid"]] = sp_count.get(o["sp_taxid"], 0) + 1
    for o in origs:
        if o["sp_taxid"] is not None and sp_count[o["sp_taxid"]] > 1 and o["anage_taxid"] is not None and o["anage_taxid"] != o["sp_taxid"]:
            review.append({"review_type": "build1_taxid_conflict", "ncbi_taxid": o["anage_taxid"], "scientific_name": o["s"].species,
                           "refseq": o["s"].refseq, "source": "MitoAge_build1", "source_name": o["s"].species,
                           "detail": f"Build 1 record taxid {o['rec_taxid']} (species {o['sp_taxid']}) is shared with another build 1 row; the AnAge name match gives species taxid {o['anage_taxid']}, used as this row's species."})
            o["eff"] = o["anage_taxid"]
    seen = set()
    for o in origs:
        if o["eff"] is None: raise ValueError(f"{o['s'].species}: no NCBI species for the cited record")
        if o["eff"] in seen: raise ValueError(f"two build 1 species resolve to NCBI species {o['eff']}")
        seen.add(o["eff"])
    orig_taxa = set()
    for o in origs:
        r = record_by_acc.get(o["s"].refseq.split(".")[0])
        for t in (o["sp_taxid"], o["eff"], o["anage_taxid"], r.species_taxid if r else None):
            if t is not None: orig_taxa.add(t)
    orig_hagrids = {o["entry"].hagrid for o in origs if o["entry"]}
    orig_species = {o["sp_taxid"] for o in origs}
    for t in [t for t in choice if t in orig_taxa and t not in orig_species]:
        o = next((x for x in origs if x["anage_taxid"] == t or x["eff"] == t), None)
        review.append({"review_type": "refseq_species_excluded_via_build1_link", "ncbi_taxid": t,
                       "scientific_name": (tax.lineage(t) or {}).get("name"), "refseq": choice[t][0].header.version, "source": None, "source_name": None,
                       "detail": f"RefSeq species linked to build 1 row {o['id'] if o else 'undefined'} ({o['s'].species if o else 'undefined'}, {o['s'].refseq if o else 'undefined'}) through its AnAge entry, while that row's cited record resolves to NCBI species {o['sp_taxid'] if o else 'undefined'}; not added as a new species."})

    taken = {s.refseq for s in st.raw.kept}
    anage_obs_by_hagrid = {o.source_record_id: o for o in st.obs if o.source == "AnAge"}
    for o in origs:
        s, cited = o["s"], o["s"].refseq
        acc = cited.split(".")[0]
        version = change = retired = None
        if cited in by_version: version = cited
        elif acc in by_accession:
            version = by_accession[acc].version; change = (cited, version, "upstream_version_change")
        else:
            succ = successor_for(choice.get(o["eff"]), taken)
            if succ:
                version = succ.header.version; taken.add(version); change = (cited, version, "refseq_successor")
            else:
                retired = (f"RefSeq {cited} is not in {RS} (suppressed) and the release has no other complete mitogenome of NCBI species "
                           f"{o['eff']} ({(tax.lineage(o['eff']) or {}).get('name', '?')}); no longer qualifies under the RefSeq-only rule")
        lo = None
        an = anage_obs_by_hagrid.get(o["entry"].hagrid) if o["entry"] else None
        group = ls.group_of_build1(s.phylum, s.cls)
        if an and an.value_years is not None: lo = an
        else:
            lo = ls.choose([x for x in by_taxid.get(o["eff"], []) if x.source != "AnAge"], group)
            if lo is None:
                why = (f"{'AnAge build 15 entry ' + o['entry'].hagrid + ' has no maximum longevity' if o['entry'] else 'no AnAge build 15 entry'}"
                       f" and no secondary lifespan source applies ({group}: {s.cls})")
                retired = f"{retired}; {why}" if retired else f"{why}; no longer qualifies (no curated lifespan)"
        o.update(version=version, change=change, retired=retired, lifespan=lo)

    cands = []
    for t, recs in choice.items():
        if t in orig_taxa: continue
        group = ls.taxon_group(recs[0].header.lineage)
        obs = by_taxid.get(t, [])
        chosen = ls.choose(obs, group)
        if chosen is None: continue
        hagrid = next((x.source_record_id for x in obs if x.source == "AnAge"), None)
        if hagrid and hagrid in orig_hagrids:
            review.append({"review_type": "anage_entry_linked_to_build1", "ncbi_taxid": t, "scientific_name": (tax.lineage(t) or {}).get("name"),
                           "refseq": recs[0].header.version, "source": "AnAge", "source_name": None,
                           "detail": f"AnAge HAGRID {hagrid} is already linked to a build 1 species (taxonomic split or reassignment); not added."})
            continue
        cands.append({"t": t, "rec": recs[0], "group": group, "chosen": chosen, "obs": obs, "hagrid": hagrid})
    return {"origs": origs, "cands": cands, "review": review, "choice": choice, "by_taxid": by_taxid, "RS": RS}


# --- stage 4: genomes, changelog, tables --------------------------------------------------------
def release_texts(st: State, versions: list[str]) -> tuple[dict, list]:
    want = set(versions)
    texts = {h.version: t for h, t in scan(st.inp.release, keep=lambda h: h.version in want) if t is not None}
    missing = sorted(want - set(texts))
    if missing: raise ValueError(f"RefSeq release: records not found: {', '.join(missing)}")
    contig = sorted(v for v, t in texts.items() if "\nORIGIN" not in t)
    by_version = {h.version: h for h in st.headers}
    for v in contig:
        p = st.inp.genbank / "contig" / f"{v}.gb"
        if not p.is_file():
            raise FileNotFoundError(f"{v} has only a CONTIG line in the release: `mitoage download contig` fetches its sequence")
        t = p.read_text()
        rec = gm.parse_genbank(t)
        if rec.accession_version != v or len(rec.sequence) != by_version[v].length:
            raise ValueError(f"efetch record for {v} is {rec.accession_version} with {len(rec.sequence)} bp (release: {by_version[v].length} bp)")
        texts[v] = t
    return texts, contig


def primary_genome(st: State, text: str, phylum: str):
    pa = primary_accession(text)
    if not pa: return None
    folder = st.inp.genbank / "primary"
    files = [p for p in folder.glob("*.gb") if p.name == f"{pa}.gb" or p.name.startswith(f"{pa}.")] if folder.is_dir() else []
    if len(files) != 1: return None
    return files[0].stem, gm.annotate(files[0].read_text(), gm.fallback_table(phylum), phylum == "Chordata")


def assemble(st: State, pl_: dict, use_primary: bool = True) -> dict:
    tax, inp, RS = st.tax, st.inp, pl_["RS"]
    date = inp.build_date
    changelog, review = [], list(pl_["review"])
    def log(c, src, ver):
        changelog.append({"build": 2, "changed_at": date, "table": c["table"], "row_key": c["row_key"], "field": c["field"],
                          "old_value": c["old_value"], "new_value": c["new_value"], "change_type": c["change_type"],
                          "diff_class": c["diff_class"], "reason": c["reason"], "upstream_source": src, "upstream_version": ver, "scope": c["scope"]})
    def sp(sid, fld, old, new, c, src, ver):
        log({"table": "species", "row_key": str(sid), "field": fld, "old_value": None if old is None else ls.js_str(old),
             "new_value": None if new is None else ls.js_str(new), **c}, src, ver)

    origs, cands = pl_["origs"], pl_["cands"]
    wanted = [o["version"] for o in origs if o["version"] and not o["retired"]] + [c["rec"].header.version for c in cands]
    st.step(f"genomes: {len(wanted)} records from the release")
    texts, contig = release_texts(st, wanted)
    anage_by_hagrid = {e.hagrid: e for e in st.anage}
    species, stats, usage, features, observations = [], [], [], {}, []
    retired_rows, status_of, kinds, flags, extra = [], {}, {}, {}, {}
    kind = lambda i, k: kinds.setdefault(i, set()).add(k)
    flag = lambda i, f: flags.setdefault(i, []).append(f)
    primary_used = 0
    raw = st.raw
    patched: dict[int, list[str]] = {}       # live species -> protein genes an annotation patch added

    st.step(f"build 1 species: status and differences ({len(origs)})")
    for o in origs:
        s, sid, b1 = o["s"], o["id"], st.base[o["id"]]
        live = not o["retired"]
        g, g_text = o["g1"], st.cited_text[sid]
        if live:
            g_text = texts[o["version"]]
            g = gm.annotate(g_text, gm.fallback_table(s.phylum), s.phylum == "Chordata")
            if g.record.accession_version != o["version"]: raise ValueError(f"release record {g.record.accession_version} != {o['version']}")
            rows = gm.genome_rows(sid, g)
            patched[sid] = g.patched
        else:
            bs, bu = build1_rows(st, sid)
            rows = {"stats": bs, "usage": bu, "features": [dict(f) for f in b1["features"]], "non_coding": b1["non_coding"] if b1["non_coding"] is not None else 0}
        primary = None
        if not o["change"] and use_primary:
            primary = primary_genome(st, g_text, s.phylum)
            if primary: primary_used += 1
        csv_comp = {r: m[sid] for r, m in raw.composition.items() if m.get(sid)}
        csv_cod = {r: m[sid] for r, m in raw.codons.items() if sid in m}
        d = dmod.diff_genome({"id": sid, "csv_composition": csv_comp, "csv_codons": csv_cod, "emitted_stats": rows["stats"],
                              "emitted_codons": rows["usage"], "evidence": g, "evidence_text": g_text, "primary": primary,
                              "acc": o["change"], "baseline_features": b1["features"], "emitted_features": rows["features"]})
        rec_version = o["version"] if live else s.refseq
        for c in d["changes"]:
            log(c, RS, rec_version)
            if c["table"] != "species": kind(sid, "genome" if o["change"] else "correction" if c["change_type"] == "correction" else "annotation")
        for k in d["removed"]:
            table, _, region = k.split(":", 2)[0], None, k.split(":")[2]
            retired_rows.append(k)
            if table == "region_stats":
                rows["stats"].append(gm.make_stats(sid, region, b1mod.fixed_counts(csv_comp[region][0]), "mitoage"))
            elif table == "codon_usage":
                rows["usage"].append(gm.usage_row(sid, region, csv_cod[region], g.transl_table, "mitoage"))
            elif table == "genome_features":
                fk = k[len(f"genome_features:{sid}:"):]
                keys = gm.feature_keys(b1["features"])
                if fk in keys: rows["features"].append(dict(b1["features"][keys.index(fk)]))
        rows["features"].sort(key=lambda f: (f["start"], f["end"]))
        if live and d["removed"]: flag(sid, f"rows_not_in_current_record:{len(d['removed'])}")

        csv_total = raw.composition["total_mtDNA"][sid][0]
        total = next(r for r in rows["stats"] if r["region"] == "total_mtDNA")
        ch = o["change"]
        acc_class = None
        if ch:
            acc_class = dmod.cls("value_update", ch[2],
                                 f"RefSeq {ch[0]} was suppressed; switched to the live RefSeq mitogenome of the same species {ch[1]} (species-level record, newest full-length)"
                                 if ch[2] == "refseq_successor" else f"RefSeq record replaced by a new sequence version {ch[0]} -> {ch[1]}", "row")
            sp(sid, "refseq", ch[0], ch[1], acc_class, RS, ch[1]); kind(sid, "genome")
        if live and csv_total.length != total["length_bp"]:
            sp(sid, "mtdna_size_bp", csv_total.length, total["length_bp"], acc_class or dmod.cls(
                "correction", "build1_disagrees_with_cited_record", "Genome length differs from the cited record", "row"), RS, rec_version)
        gc_now = gm.js_round(total["gc_percent"], 3) if live else b1["gc_percent"]
        if live and acc_class and gc_now != b1["gc_percent"]: sp(sid, "gc_percent", b1["gc_percent"], gc_now, acc_class, RS, rec_version)
        nonc = rows["non_coding"] if live else b1["non_coding"]
        if csv_total.non_coding != nonc:
            if acc_class: c = acc_class
            elif d["aggregates"]:
                c = dmod.cls("value_update" if all(a[1]["change_type"] == "value_update" for a in d["aggregates"]) else "correction",
                             "derived_from_region_differences",
                             f"Non-coding length = genome minus CDS/rRNA/tRNA coverage; differs because aggregate regions differ ({', '.join(a[0] for a in d['aggregates'])})",
                             "rule" if all(a[1]["scope"] == "rule" for a in d["aggregates"]) else "row")
            else:
                c = dmod.cls("correction", "build1_disagrees_with_cited_record",
                             "Build-1 non-coding length not reproducible from the cited record (genome minus union of CDS, rRNA, tRNA); no aggregate-region difference explains it", "row")
            sp(sid, "non_coding_bp", csv_total.non_coding, nonc, c, RS, rec_version)
            kind(sid, "genome" if ch else "correction" if c["change_type"] == "correction" else "annotation")
        rec_taxid_now = (None if g.taxid is None else tax.current(g.taxid)) if live else o["rec_taxid"]
        if o["g1"].taxid != rec_taxid_now:
            sp(sid, "ncbi_taxid", o["g1"].taxid, rec_taxid_now, acc_class or dmod.cls(
                "value_update", "ncbi_taxid_merged", f"NCBI merged taxid {o['g1'].taxid} into {rec_taxid_now}", "row"), "NCBI Taxonomy", st.inp.taxonomy_label)
            kind(sid, "ncbi_taxid")
        lo = o["lifespan"]; e = o["entry"]
        max_lon = lo.value_years if lo else s.max_longevity
        if lo and max_lon != s.max_longevity:
            u = lambda v: "undefined" if e is None else v
            sp(sid, "max_longevity_years", s.max_longevity, max_lon, dmod.cls(
                "value_update", "anage_longevity_revised" if lo.source == "AnAge" else "secondary_source_fills_missing_anage",
                f"AnAge build 15 revised maximum longevity (quality {u(e and e.data_quality)}, sample size {u(e and e.sample_size)}, origin {u(e and e.specimen_origin)})"
                if lo.source == "AnAge" else f"AnAge build 15 has no value; {lo.source} fills it under the per-taxon priority", "row"), lo.source, lo.source_version)
            kind(sid, "lifespan")
        ncbi_name = (tax.lineage(o["eff"]) or {}).get("name")
        rec_sp_name = None if o["sp_taxid"] is None else (tax.lineage(o["sp_taxid"]) or {}).get("name")
        name, previous = s.species, []
        nc = name_change(s.species, ncbi_name)
        if nc in RENAME_CLASSES:
            syn = s.species in (tax.names.get(o["eff"]) or {}).get("synonyms", [])
            sp(sid, "scientific_name", s.species, ncbi_name, dmod.cls(
                "renamed", f"ncbi_{nc}", f"{nc}; build-1 name {'is' if syn else 'is not'} listed as an NCBI synonym; AnAge build 15 name: {e.scientific_name if e else '-'}", "row"),
               "NCBI Taxonomy", st.inp.taxonomy_label)
            name = ncbi_name; previous.append(s.species); kind(sid, "name")
        rec_nc = name_change(s.species, rec_sp_name)
        if rec_nc == "different_species_name" or nc == "different_species_name":
            flag(sid, f"cited_record_resolves_to_other_species:{rec_sp_name}")
        if e and name_change(name, e.scientific_name) == "different_species_name": flag(sid, f"anage_name_differs:{e.scientific_name}")
        if o["eff"] != o["sp_taxid"]: flag(sid, f"ncbi_taxid_shared_in_build1:{o['sp_taxid']}")
        genus = name.split(" ")[0]
        common = (e.common_name if e and e.common_name else s.common_name)
        taxo = {"kingdom": s.kingdom, "phylum": (e.phylum if e and e.phylum else s.phylum), "class": s.cls,
                "order": (e.order if e and e.order else s.order), "family": (e.family if e and e.family else s.family), "genus": genus}
        if genus != s.genus:
            sp(sid, "genus", s.genus, genus, dmod.cls("renamed", "follows_scientific_name", "Genus follows the renamed scientific name", "row"),
               "NCBI Taxonomy", st.inp.taxonomy_label)
        for fld, old, now in (("common_name", s.common_name, common), ("phylum", s.phylum, taxo["phylum"]),
                              ("order", s.order, taxo["order"]), ("family", s.family, taxo["family"])):
            if old == now: continue
            sp(sid, fld, old, now, dmod.cls("value_update", "anage_taxonomy_revision", "AnAge build 15 taxonomy/common-name revision", "row"),
               "AnAge", pins.ANAGE_VERSION)
            kind(sid, "common_name" if fld == "common_name" else "taxonomy")
        if o["retired"]:
            no_record = not o["version"]
            sp(sid, "record_status", "original", "retired", dmod.cls(
                "retired", "refseq_suppressed_no_successor" if no_record else "no_curated_lifespan", o["retired"], "row"),
               RS if no_record else "AnAge", s.refseq if no_record else pins.ANAGE_VERSION)
        lifespan = (lifespan_of(lo) if lo else {"source": "MitoAge_build1", "source_version": "MitoAge build 1 (2015-07-01)",
                    "source_record_id": str(sid), "quality": None, "sample_size": None, "record_type": "compiled", "license": pins.BUILD1_LICENSE})
        if lo is None: flag(sid, "lifespan_from_build1_only")
        species.append({"id": sid, "slug": b1["slug"], "refseq": o["version"] if live else s.refseq, "scientific_name": name,
                        "common_name": common, **taxo, "max_longevity_years": max_lon,
                        "mtdna_size_bp": total["length_bp"] if live else b1["mtdna_size"], "non_coding_bp": nonc, "gc_percent": gc_now,
                        "legacy_id": b1["legacy_id"], "traits": b1["traits"], "ncbi_taxid": rec_taxid_now, "ncbi_species_taxid": o["eff"],
                        "lifespan": lifespan, "previous": previous, "retired_reason": o["retired"]})
        status_of[sid] = "retired" if o["retired"] else "renamed" if previous else "original"
        stats += rows["stats"]; usage += rows["usage"]; features[sid] = rows["features"]
        e0, how = st.anage_match[sid]
        extra[sid] = {"anage_hagrid": e0.hagrid if e0 else None, "anage_match": how, "anage_max": e0.max_longevity if e0 else None,
                      "transl_table": g.transl_table}
        an = None
        if e:
            an = next((x for x in reversed(st.obs) if x.source == "AnAge" and x.source_record_id == e.hagrid), None)
        here = []
        if an and an.value_years is not None: here.append(ls.replace(an, match_method="build1_hagrid_link"))
        here += [x for x in pl_["by_taxid"].get(o["eff"], []) if x.source != "AnAge" and x.value_years is not None]
        for x in here:
            observations.append(obs_row(sid, o["eff"], x, x is lo or (lo is not None and x.source_record_id == lo.source_record_id and x.source == lo.source)))
        observations.append({"species_id": sid, "ncbi_species_taxid": o["eff"], "source": "MitoAge_build1", "source_version": "MitoAge build 1 (2015-07-01)",
                             "source_record_id": str(sid), "source_name": s.species, "value_years": s.max_longevity, "record_type": "compiled",
                             "quality": None, "sample_size_text": None, "sex": None, "detail": "build-1 max_longevity_years (AnAge build 12-era)",
                             "reference": None, "license": pins.BUILD1_LICENSE, "match_method": "build1_row", "is_chosen": lo is None})

    st.step(f"new species ({len(cands)})")
    used_slugs = {s["slug"] for s in species}; used_ids = {s["id"] for s in species}
    for c in sorted(cands, key=lambda x: x["t"]):
        sid = mint_id(c["t"])
        if sid in used_ids: raise ValueError(f"minted id {sid} collides")
        used_ids.add(sid)
        lin = tax.lineage(c["t"])
        if not lin: raise ValueError(f"no NCBI lineage for species {c['t']}")
        ver = c["rec"].header.version
        g = gm.annotate(texts[ver], gm.fallback_table(lin["phylum"]), lin["phylum"] == "Chordata")
        rows = gm.genome_rows(sid, g)
        patched[sid] = g.patched
        entry = anage_by_hagrid.get(c["hagrid"]) if c["hagrid"] else None
        name = lin["name"]
        slug = slugify(name)
        if slug in used_slugs: slug = f"{slug}-{slugify(ver)}"
        used_slugs.add(slug)
        total = next(r for r in rows["stats"] if r["region"] == "total_mtDNA")
        traits = ({**entry.traits, "source": f"AnAge build 15 (HAGRID {entry.hagrid})"} if entry else {**{k: None for k in ls.TRAITS}, "source": None})
        common = (entry.common_name if entry and entry.common_name else None) or c["chosen"].source_common_name or \
                 ((tax.names.get(c["t"]) or {}).get("genbank_common")) or ""
        hc = ls.harmonised_class(c["rec"].header.lineage, lin["class"])
        species.append({"id": sid, "slug": slug, "refseq": ver, "scientific_name": name, "common_name": common, "kingdom": "Animalia",
                        "phylum": lin["phylum"], "class": hc if hc is not None else (entry.cls if entry else ""), "order": lin["order"],
                        "family": lin["family"], "genus": lin["genus"] or name.split(" ")[0], "max_longevity_years": c["chosen"].value_years,
                        "mtdna_size_bp": total["length_bp"], "non_coding_bp": rows["non_coding"], "gc_percent": gm.js_round(total["gc_percent"], 3),
                        "legacy_id": None, "traits": traits, "ncbi_taxid": c["rec"].taxid_current, "ncbi_species_taxid": c["t"],
                        "lifespan": lifespan_of(c["chosen"]), "previous": [], "retired_reason": None})
        status_of[sid] = "new"
        stats += rows["stats"]; usage += rows["usage"]; features[sid] = rows["features"]
        extra[sid] = {"anage_hagrid": entry.hagrid if entry else None, "anage_match": "ncbi name match" if entry else "none",
                      "anage_max": entry.max_longevity if entry else None, "transl_table": g.transl_table}
        h = c["rec"].header
        added = dmod.cls("added", "new_species",
                         f"Complete RefSeq mitogenome ({ver}, {c['rec'].completeness}{', ' + str(h.length) + ' bp' if h.length else ''}) and a maximum lifespan from "
                         f"{c['chosen'].source} ({c['group']} priority: {'AnAge' if c['chosen'].source == 'AnAge' else 'AnAge has no value'})", "row")
        log({"table": "species", "row_key": str(sid), "field": "(row)", "old_value": None, "new_value": f"{name} {ver}", **added}, RS, ver)
        for table, n in (("region_stats", len(rows["stats"])), ("codon_usage", len(rows["usage"])), ("genome_features", len(rows["features"]))):
            log({"table": table, "row_key": f"{sid}:*", "field": "(rows)", "old_value": None, "new_value": f"{n} rows", **added,
                 "reason": f"Rows of new species {name}"}, RS, ver)
        for x in c["obs"]:
            if x.value_years is not None: observations.append(obs_row(sid, c["t"], x, x is c["chosen"]))

    # protein genes a record leaves out (gene_annotation.py, the site's D17): a patched gene is flagged and logged;
    # a missing gene carries its reason as a flag, and one with no recorded reason is flagged as unexplained and
    # listed for review. An entry of the lists that no live species uses is refused, as the site's code refuses it.
    from .gene_annotation import ANNOTATION_PATCHES, GENE_GAPS, LINEAGE_GENE_GAPS
    live = [s for s in species if status_of.get(s["id"]) != "retired" and s["id"] in patched]
    live_records = {s["refseq"] for s in live}
    for p in ANNOTATION_PATCHES:
        if p.accession not in live_records: raise ValueError(f"annotation patch {p.accession} {p.gene}: no live species uses that record")
    gaps_used = set()
    for s in live:
        sid = s["id"]
        for gene in patched[sid]:
            p = next(x for x in ANNOTATION_PATCHES if x.accession == s["refseq"] and x.gene == gene)
            flag(sid, f"gene_annotated_by_mitoage:{gene}"); kind(sid, "annotation")
            log({"table": "genome_features", "row_key": f"{sid}:{gene}#1", "field": "(row)", "old_value": None,
                 "new_value": f"{gene} CDS {p.location}", "change_type": "correction", "diff_class": "annotation_patch",
                 "reason": p.evidence, "scope": "row"}, RS, s["refseq"])
        present = {f["name"] for f in features.get(sid, []) if f["type"] == "CDS"}
        for gene in gm.PROTEIN_GENES:
            own = next((x for x in GENE_GAPS if x.accession == s["refseq"] and x.gene == gene), None)
            if gene in present:
                if own: raise ValueError(f"gene gap {own.accession} {gene}: the build now has {gene} for {s['scientific_name']}")
                continue
            gap = own or next((x for x in LINEAGE_GENE_GAPS if x.gene == gene and s["phylum" if x.rank == "phylum" else "class"] == x.taxon), None)
            if own: gaps_used.add(own)
            if gap:
                flag(sid, f"gene_{gap.reason}:{gene}")
            else:
                flag(sid, f"gene_unexplained_absence:{gene}")
                review.append({"review_type": "protein_gene_not_annotated", "ncbi_taxid": s["ncbi_species_taxid"],
                               "scientific_name": s["scientific_name"], "refseq": s["refseq"], "source": None, "source_name": None,
                               "detail": f"{gene} is not annotated in the record and no reason is recorded in pipeline/sources/gene-annotation.ts"})
    for x in GENE_GAPS:
        if x not in gaps_used: raise ValueError(f"gene gap {x.accession} {x.gene}: no live species uses that record without the gene")

    # review flags, possible duplicates, the review list
    by_name = {s["scientific_name"]: s for s in species}
    for a, b, why in POSSIBLE_DUPLICATES:
        sa, sb = by_name.get(a), by_name.get(b)
        if not sa or not sb: continue
        flag(sa["id"], f"possible_duplicate_species:{b}"); flag(sb["id"], f"possible_duplicate_species:{a}")
        review.append({"review_type": "possible_duplicate_species", "ncbi_taxid": None, "scientific_name": f"{a} ({sa['id']}); {b} ({sb['id']})",
                       "refseq": f"{sa['refseq']}; {sb['refseq']}", "source": None, "source_name": None, "detail": f"{why}; both kept and flagged (owner default h)"})
    choice = pl_["choice"]
    genera = {((tax.lineage(t) or {}).get("name") or "").split(" ")[0].lower() for t in choice}
    in_build = {s["ncbi_species_taxid"] for s in species}
    for m in sorted(st.matches.values(), key=lambda m: (m.source + m.source_name, m.speccode or 0)):
        genus = m.source_name.lower().split(" ")[0]
        if m.status == "unmatched" and genus in genera:
            review.append({"review_type": "source_name_unmatched", "ncbi_taxid": None, "scientific_name": None, "refseq": None, "source": m.source,
                           "source_name": m.source_name, "detail": f"no NCBI scientific name or synonym; the release has complete mitogenomes in genus {m.source_name.split(' ')[0]}"})
        elif m.status in ("conflict_secondary", "ambiguous", "ambiguous_many_to_one"):
            ids = m.candidates + ([m.species_taxid] if m.species_taxid is not None else [])
            if not any(t in choice for t in ids): continue
            review.append({"review_type": f"source_name_{m.status}", "ncbi_taxid": m.species_taxid,
                           "scientific_name": (tax.lineage(m.species_taxid) or {}).get("name") if m.species_taxid else None, "refseq": None,
                           "source": m.source, "source_name": m.source_name,
                           "detail": f"match_method={m.method}; candidate taxids={';'.join(map(str, m.candidates)) or '-'}; not used"})
            for t in ids:
                s = next((x for x in species if x["ncbi_species_taxid"] == t), None)
                if s and t in in_build: flag(s["id"], f"source_match_in_review:{m.source}:{m.source_name}")

    # rules logged once
    def rule(table, fld, old, new, dc, reason, src, ver, key="*"):
        changelog.append({"build": 2, "changed_at": date, "table": table, "row_key": key, "field": fld, "old_value": old, "new_value": new,
                          "change_type": "correction", "diff_class": dc, "reason": reason, "upstream_source": src, "upstream_version": ver, "scope": "rule"})
    for f, col in (("a", "T"), ("c", "G"), ("g", "C"), ("t", "A")):
        rule("region_stats", f, f'CSV column "{col}"', f'true RefSeq-strand {f.upper()} count (= CSV "{col}")', "rule_label_complement",
             "Build-1 CSV base columns are complement-labelled (G<->C, A<->T); verified on total_mtDNA for every species; GC% unaffected, GC and AT skews flip sign",
             "MitoAge build-1 CSV export vs RefSeq sequences", "pipeline rule (fixCsvLabels)")
    rule("region_stats", "length_bp", "start..end span of segmented features", "sum of annotated segments only", "rule_full_span",
         "Build 1 counted the full start..end span of join() features (e.g. avian/turtle ND3 frameshift nucleotide); build 2 counts annotated segments only (owner default f); per-row instances carry diff_class=rule_full_span",
         "RefSeq feature tables", "pipeline rule (coverage of segments)")
    rule("species", "record_status", "in CSV export (923rd row)", "excluded", "rule_exclusion",
         "Loxioides bailleui NC_025626.1 is 9,243 bp (not a complete mitogenome); not in the published 922; excluded from every build",
         "NCBI nuccore", "NC_025626.1", "HAGRID 1598")
    changelog.append({"build": 2, "changed_at": date, "table": "lifespan_observations", "row_key": "*", "field": "(table)", "old_value": None,
                      "new_value": f"{len(observations)} rows", "change_type": "added", "diff_class": "new_table",
                      "reason": "Every lifespan value from every source kept in long format (decisions section 2); new in build 2",
                      "upstream_source": "AnAge, FishBase, SeaLifeBase, USGS BBL", "upstream_version": "see source_version", "scope": "rule"})

    # provenance
    touched: dict[str, bool] = {}
    for c in changelog:
        if c["row_key"] == "*" or c["row_key"].endswith(":*"): continue
        k = f"{c['table']}:{c['row_key']}"
        touched[k] = touched.get(k, False) or (c["field"] == "(row)" and c["old_value"] is None)
    retired_keys = set(retired_rows)
    modified = {int(k.split(":")[1]) for k in touched if k.split(":")[1].isdigit()}   # "species:HAGRID 1598" names no row
    def prov(added, status, mod, ver):
        return {"added_in_build": added, "added_at": pins.BUILD1_DATE if added == 1 else date, "last_modified_build": 2 if mod else None,
                "last_modified_at": date if mod else None, "record_status": status, "source_version": ver}
    for s in species:
        st0 = status_of[s["id"]]
        mod = s["id"] in modified
        stt = "updated" if st0 == "original" and mod else st0
        s["review_flags"] = sorted(set(flags.get(s["id"], [])))
        s["change_kinds"] = sorted(kinds.get(s["id"], set()))
        s.update(prov(2 if stt == "new" else 1, stt, False if stt == "new" else mod, f"{s['refseq']}; {s['lifespan']['source_version']}"))
    sp_status = {s["id"]: s["record_status"] for s in species}
    def child(table, sid, key, ver):
        sst = sp_status[sid]; k = f"{table}:{key}"; t = touched.get(k)
        if sst == "new": return prov(2, "new", False, ver)
        if sst == "retired" or k in retired_keys: return prov(2 if t else 1, "retired", t is not None and not t, ver)
        if t: return prov(2, "new", False, ver)
        return prov(1, "updated" if t is not None else "original", t is not None, ver)
    ref = {s["id"]: s["refseq"] for s in species}
    for r in stats: r.update(child("region_stats", r["species_id"], f"{r['species_id']}:{r['region']}", f"{ref[r['species_id']]} (build 1 CSV)" if r["source"] == "mitoage" else ref[r["species_id"]]))
    for u in usage: u.update(child("codon_usage", u["species_id"], f"{u['species_id']}:{u['region']}", f"{ref[u['species_id']]} (build 1 CSV)" if u["source"] == "mitoage" else ref[u["species_id"]]))
    feat_rows = []
    for sid in sorted(features):
        fs = features[sid]
        for k, f in zip(gm.feature_keys(fs), fs):
            feat_rows.append({"species_id": sid, "feature_key": k, **f, **child("genome_features", sid, f"{sid}:{k}", ref[sid])})
    for o in observations: o.update(prov(2, "new", False, o["source_version"]))

    # deterministic order
    species.sort(key=lambda s: s["id"])
    ro = {r: i for i, r in enumerate(gm.REGIONS)}; co = {r: i for i, r in enumerate(gm.CODON_REGIONS)}
    stats.sort(key=lambda r: (r["species_id"], ro[r["region"]]))
    usage.sort(key=lambda r: (r["species_id"], co[r["region"]]))
    observations.sort(key=lambda o: (o["species_id"], SRC_ORDER.index(o["source"]), o["source_record_id"]))
    key_num = lambda k: -1 if k == "*" or k.startswith("HAGRID") else int(k.split(":")[0])
    changelog.sort(key=lambda c: (key_num(c["row_key"]), c["table"], c["row_key"], c["field"], c["diff_class"], c["old_value"] or ""))
    review.sort(key=lambda r: (r["review_type"], r["source"] or "", r["source_name"] or r["scientific_name"] or ""))
    return {"species": species, "stats": stats, "usage": usage, "features": feat_rows, "observations": observations,
            "changelog": changelog, "review": review, "extra": extra,
            "report": {"primary_records_used": primary_used, "candidates": len(cands), "retired_child_rows": len(retired_rows),
                       "contig_only_fetched": contig}}


def lifespan_of(o) -> dict:
    return {"source": o.source, "source_version": o.source_version, "source_record_id": o.source_record_id, "quality": o.quality,
            "sample_size": o.sample_size_text if o.source == "AnAge" else None, "record_type": o.record_type, "license": o.license}


def obs_row(sid: int, taxid: int, o, chosen: bool) -> dict:
    return {"species_id": sid, "ncbi_species_taxid": taxid, "source": o.source, "source_version": o.source_version,
            "source_record_id": o.source_record_id, "source_name": o.source_name, "value_years": o.value_years, "record_type": o.record_type,
            "quality": o.quality, "sample_size_text": o.sample_size_text, "sex": o.sex, "detail": o.detail, "reference": o.reference,
            "license": o.license, "match_method": o.match_method, "is_chosen": chosen}
