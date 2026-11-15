"""Lifespan sources, their matching to NCBI species and the per-taxon precedence of build 2 (site
repository: pipeline/sources/{anage,fishbase,usgs-bbl,euring,observation}.ts, pipeline/steps/lifespan.ts,
plans/DATA_BUILD2_DECISIONS.md section 2).

Values read: AnAge "Maximum longevity (yrs)" (record type from "Specimen origin": wild, captivity ->
captive, otherwise observed; quality from "Data quality"); FishBase and SeaLifeBase popchar.tmax
(observed; tmaxQuality -1 -> questionable, else acceptable), species.LongevityWild (wild, acceptable),
species.LongevityCaptive (captive, acceptable), popgrowth.tmax (estimated, low); USGS BBL "Minimum age at
encounter" (wild, acceptable; the oldest banded bird, a minimum). EURING: the list of 2023-09-01 (euring.py), each
record wild and acceptable, a minimum age.

Matching order: NCBI scientific name, NCBI synonym, the source's own synonym table (FishBase, SeaLifeBase;
misapplied names never), binomial of a trinomial. A name that gives several species is ambiguous; several
records of one source on one species keep the one with the strictly best method, and a tie makes all of
them ambiguous. Ambiguous records are not used.

Precedence: AnAge always first. Without an AnAge value: fish take FishBase, aquatic invertebrates
SeaLifeBase (only species flagged fresh, brackish or salt water), birds EURING then USGS BBL, other
vertebrates nothing. Within a source: reported values that are not questionable, then reported
questionable ones, then estimates; the maximum within the first tier that has values; at equal maxima a
value with a stated specimen origin (wild, captive) beats one without (observed), then the first listed.
"""
import csv
import io
import math
import re
import zipfile
from dataclasses import dataclass, replace
from pathlib import Path

import polars as pl

from . import pins
from .taxonomy import name_key

PRECEDENCE = {"fish": ["AnAge", "FishBase"], "invertebrate": ["AnAge", "SeaLifeBase"],
              "bird": ["AnAge", "EURING", "USGS_BBL"], "other_vertebrate": ["AnAge"]}
METHOD_RANK = {"ncbi_scientific_name": 0, "ncbi_synonym": 1, "fishbase_synonym": 2, "sealifebase_synonym": 2,
               "binomial_of_trinomial": 3}
FISH_CLASSES = {"Actinopterygii", "Chondrichthyes", "Myxini", "Sarcopterygii", "Cephalaspidomorphi"}


@dataclass
class Obs:
    source: str
    source_version: str
    source_record_id: str
    source_name: str
    source_common_name: str | None
    value_years: float | None
    record_type: str
    quality: str | None
    sample_size_text: str | None
    sex: str | None
    detail: str | None
    reference: str | None
    license: str
    speccode: int | None = None
    aquatic: bool | None = None
    species_taxid: int | None = None     # set once matched
    match_method: str | None = None


def collapse(s) -> str:
    return re.sub(r"\s+", " ", "" if s is None else str(s)).strip()


def js_str(v) -> str:
    """String(v) as JavaScript writes a number read from parquet (12.0 -> "12")."""
    if isinstance(v, bool): return "true" if v else "false"
    if isinstance(v, float):
        if v.is_integer() and abs(v) < 1e21: return str(int(v))
        r = repr(v)
        m = re.match(r"^(-?[\d.]+)e([+-])0*(\d+)$", r)
        return f"{m.group(1)}e{m.group(2)}{m.group(3)}" if m else r
    return str(v)


def num(v) -> float | None:
    if v is None: return None
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return float(v) if math.isfinite(v) else None
    s = str(v).strip()
    if s == "": return None
    try:
        x = float(s)
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def text(v) -> str | None:
    return None if v is None else js_str(v)


def quality(s) -> str | None:
    v = (s or "").strip().lower()
    return v if v in ("high", "acceptable", "questionable", "low") else None


def read_tsv(text_: str, delimiter: str) -> list[dict]:
    rows = [r for r in csv.reader(io.StringIO(text_), delimiter=delimiter) if len(r) > 1 or (r and r[0] != "")]
    header = [h.strip() for h in rows[0]]
    return [{h: (r[i] if i < len(r) else "").strip() for i, h in enumerate(header)} for r in rows[1:]]


# --- AnAge ------------------------------------------------------------------------------------
@dataclass
class AnAgeEntry:
    hagrid: str
    scientific_name: str
    common_name: str
    cls: str
    kingdom: str
    phylum: str
    order: str
    family: str
    max_longevity: float | None
    data_quality: str
    sample_size: str
    specimen_origin: str
    source_text: str
    references: str
    traits: dict


TRAITS = {"adult_weight_g": "Adult weight (g)", "birth_weight_g": "Birth weight (g)", "metabolic_rate_w": "Metabolic rate (W)",
          "body_temperature_k": "Temperature (K)", "female_maturity_days": "Female maturity (days)",
          "male_maturity_days": "Male maturity (days)", "gestation_days": "Gestation/Incubation (days)",
          "litter_size": "Litter/Clutch size", "litters_per_year": "Litters/Clutches per year",
          "growth_rate_per_day": "Growth rate (1/days)"}


def load_anage(path: Path) -> list[AnAgeEntry]:
    if not path.is_file():
        raise FileNotFoundError(f"no {path}: run `mitoage download lifespan` first")
    with zipfile.ZipFile(path) as z:
        name = next((n for n in z.namelist() if n.endswith("anage_data.txt")), None)
        if name is None: raise ValueError(f"{path} holds no anage_data.txt")
        rows = read_tsv(z.read(name).decode("utf-8"), "\t")
    return [AnAgeEntry(r["HAGRID"], collapse(f"{r.get('Genus', '')} {r.get('Species', '')}"), r.get("Common name", ""),
                       r.get("Class", ""), r.get("Kingdom", ""), r.get("Phylum", ""), r.get("Order", ""), r.get("Family", ""),
                       num(r.get("Maximum longevity (yrs)")), r.get("Data quality", ""), r.get("Sample size", ""),
                       r.get("Specimen origin", ""), r.get("Source", ""), r.get("References", ""),
                       {k: num(r.get(c)) for k, c in TRAITS.items()}) for r in rows]


def anage_observations(entries: list[AnAgeEntry]) -> list[Obs]:
    """Every animal entry (entries without a value are kept for name matching)."""
    return [Obs("AnAge", pins.ANAGE_VERSION, e.hagrid, e.scientific_name, e.common_name or None, e.max_longevity,
                "wild" if e.specimen_origin == "wild" else "captive" if e.specimen_origin == "captivity" else "observed",
                quality(e.data_quality), e.sample_size or None, None,
                f"specimen origin={e.specimen_origin}; source={e.source_text}", e.references or None, pins.ANAGE_LICENSE)
            for e in entries if e.kingdom == "Animalia"]


# --- FishBase, SeaLifeBase --------------------------------------------------------------------
def load_fishbase_like(folder: Path, db: str) -> tuple[list[Obs], dict[int, list[tuple[str, str]]]]:
    source = "FishBase" if db == "fb" else "SeaLifeBase"
    version, lic = pins.FISHBASE_LABEL[db], pins.FISHBASE_LICENSE[db]
    def table(t):
        p = folder / f"{t}.parquet"
        if not p.is_file(): raise FileNotFoundError(f"no {p}: run `mitoage download lifespan` first")
        return pl.read_parquet(p).to_dicts()
    species = table("species")
    by_code = {}
    for r in species:
        by_code[int(r["SpecCode"])] = (collapse(f"{str(r.get('Genus') or '').strip()} {str(r.get('Species') or '').strip()}"),
                                       text(r.get("FBname")), any(num(r.get(c)) == 1 for c in ("Fresh", "Brack", "Saltwater")))
    def base(code, **kw):
        name, common, aquatic = by_code[code]
        return Obs(source=source, source_version=version, source_name=name, source_common_name=common, license=lic,
                   speccode=code, aquatic=aquatic if source == "SeaLifeBase" else None, **kw)
    obs = []
    for r in table("popchar"):
        code, v = int(r["Speccode"]), num(r.get("tmax"))
        if v is None or code not in by_code: continue
        q = num(r.get("tmaxQuality")); q = 0 if q is None else q
        loc = r.get("Locality"); loc = "" if loc is None else js_str(loc)
        obs.append(base(code, source_record_id=f"popchar:{js_str(r['Autoctr'])}", value_years=v, record_type="observed",
                        quality="questionable" if q == -1 else "acceptable", sample_size_text=None, sex=text(r.get("Sex")),
                        detail=f"popchar tmax; locality={loc}; tmaxQuality={js_str(q)}", reference=text(r.get("PopCharRefNo"))))
    for col, rt, ref in (("LongevityWild", "wild", "LongevityWildRef"), ("LongevityCaptive", "captive", "LongevityCapRef")):
        for r in species:
            v = num(r.get(col))
            if v is None: continue
            code = int(r["SpecCode"])
            obs.append(base(code, source_record_id=f"species.{col}:{code}", value_years=v, record_type=rt, quality="acceptable",
                            sample_size_text=None, sex=None, detail=f"species.{col}", reference=text(r.get(ref))))
    for r in table("popgrowth"):
        code, v = int(r["SpecCode"]), num(r.get("tmax"))
        if v is None or code not in by_code: continue
        n = num(r.get("Number"))
        f = lambda k: "" if r.get(k) is None else js_str(r.get(k))
        obs.append(base(code, source_record_id=f"popgrowth:{js_str(r['AutoCtr'])}", value_years=v, record_type="estimated",
                        quality="low", sample_size_text=None if n is None else f"n={math.trunc(n)}", sex=text(r.get("Sex")),
                        detail=f"popgrowth tmax (from growth study; data={f('Data')}; K={f('K')}; Loo={f('Loo')})",
                        reference=text(r.get("tmaxRef"))))
    synonyms: dict[int, list[tuple[str, str]]] = {}
    for r in table("synonyms"):
        if str(r.get("TaxonLevel") or "").lower() not in ("species", "subspecies", "nominotypical"): continue
        synonyms.setdefault(int(r["SpecCode"]), []).append(
            (f"{'' if r.get('SynGenus') is None else r['SynGenus']} {'' if r.get('SynSpecies') is None else r['SynSpecies']}",
             "" if r.get("Status") is None else str(r["Status"])))
    return obs, synonyms


# --- USGS Bird Banding Laboratory -------------------------------------------------------------
def bbl_age(s: str) -> float | None:
    """'12 years 07 months' -> 12.5833."""
    m = re.match(r"^\s*(\d+)\s*years?\s*(\d+)\s*months?", s)
    return int(m.group(1)) + int(m.group(2)) / 12 if m else None


def load_usgs(path: Path) -> list[Obs]:
    if not path.is_file():
        raise FileNotFoundError(f"no {path}: run `mitoage download lifespan` first")
    rows = read_tsv(path.read_text(encoding="utf-8"), ",")
    return [Obs("USGS_BBL", pins.USGS_VERSION, f"usgs2025:{i}", collapse(r["Scientific Name"]), r["Species Name"] or None,
                bbl_age(r["Minimum age at encounter"]), "wild", "acceptable",
                "single individual (oldest banded bird; minimum age)", r["Sex"] or None,
                f"min age {r['Minimum age at encounter']}; banded {r['Banding date']} {r['Banding region']} ({r['Age at banding']}); "
                f"encountered {r['Encounter date']} {r['Encounter region']}; {r['How bird and band were obtained']}",
                "NABBP longevity table", pins.USGS_LICENSE) for i, r in enumerate(rows)]


# --- matching ---------------------------------------------------------------------------------
def rec_key(source: str, name: str, speccode) -> tuple:
    return (source, name, speccode)


def wanted_name_keys(obs: list[Obs], synonyms: dict) -> set[str]:
    keys = set()
    def add(k):
        keys.add(k); w = k.split(" ")
        if len(w) >= 3: keys.add(f"{w[0]} {w[1]}")
    for o in obs: add(name_key(o.source_name))
    for table in synonyms.values():
        for lst in table.values():
            for name, _ in lst: add(name_key(name))
    return keys


@dataclass
class Match:
    source: str
    source_name: str
    speccode: int | None
    method: str
    species_taxid: int | None
    candidates: list[int]
    status: str


def match_names(obs: list[Obs], tax, synonyms: dict) -> dict[tuple, Match]:
    recs = {}
    for o in obs:
        recs.setdefault(rec_key(o.source, o.source_name, o.speccode), o)
    out: dict[tuple, Match] = {}
    for key, r in recs.items():
        k = name_key(r.source_name)
        res = None
        if k in tax.sci: res = ("ncbi_scientific_name", tax.sci[k])
        elif k in tax.syn: res = ("ncbi_synonym", tax.syn[k])
        table = synonyms.get(r.source)
        if res is None and table is not None and r.speccode is not None:
            hits: dict[str, set[int]] = {}
            for name, status in table.get(r.speccode, []):
                sk = name_key(name)
                if status == "misapplied name" or sk == k: continue
                c = tax.sci.get(sk)
                if not c: continue
                hits.setdefault(status, set()).update(c)
            good = set()
            for st in ("synonym", "accepted name", "provisionally accepted name"): good |= hits.get(st, set())
            src = r.source.lower()
            if good: res = (f"{src}_synonym", sorted(good))
            elif hits.get("ambiguous synonym"): res = (f"{src}_ambiguous_synonym", sorted(hits["ambiguous synonym"]))
        if res is None:
            w = k.split(" ")
            if len(w) >= 3 and f"{w[0]} {w[1]}" in tax.sci: res = ("binomial_of_trinomial", tax.sci[f"{w[0]} {w[1]}"])
        if res is None:
            out[key] = Match(r.source, r.source_name, r.speccode, "unmatched", None, [], "unmatched"); continue
        ok = len(res[1]) == 1 and "ambiguous" not in res[0]
        out[key] = Match(r.source, r.source_name, r.speccode, res[0], res[1][0] if ok else None, list(res[1]), "ok" if ok else "ambiguous")
    groups: dict[tuple, list[Match]] = {}
    for m in out.values():
        if m.status == "ok": groups.setdefault((m.source, m.species_taxid), []).append(m)
    for g in groups.values():
        if len(g) <= 1: continue
        best = min(METHOD_RANK.get(m.method, 9) for m in g)
        top = [m for m in g if METHOD_RANK.get(m.method, 9) == best]
        for m in g:
            m.status = ("ok_conflict_primary" if m is top[0] else "conflict_secondary") if len(top) == 1 else "ambiguous_many_to_one"
    return out


def usable(m: Match | None) -> bool:
    return m is not None and m.status in ("ok", "ok_conflict_primary")


# --- precedence -------------------------------------------------------------------------------
def taxon_group(lineage: str) -> str:
    parts = re.split(r";\s*", lineage)
    if "Aves" in parts: return "bird"
    if re.search(r"\b(Actinopterygii|Chondrichthyes|Cyclostomata|Dipnomorpha|Actinistia)\b", lineage): return "fish"
    if "Vertebrata" not in parts: return "invertebrate"
    return "other_vertebrate"


def group_of_build1(phylum: str, cls: str) -> str:
    if cls == "Aves": return "bird"
    if cls in FISH_CLASSES: return "fish"
    if phylum != "Chordata": return "invertebrate"
    return "other_vertebrate"


def harmonised_class(lineage: str, ncbi_class: str) -> str | None:
    """The class in build 1's scheme (AnAge build 12-era labels), so scopes stay comparable with 2016."""
    p = re.split(r";\s*", lineage)
    if "Actinopterygii" in p: return "Actinopterygii"
    if any(x in p for x in ("Lepidosauria", "Testudines", "Crocodylia")): return "Reptilia"
    for c in ("Aves", "Mammalia", "Amphibia", "Chondrichthyes"):
        if c in p: return c
    if "Dipnomorpha" in p or "Actinistia" in p: return "Sarcopterygii"
    if "Hyperoartia" in p: return "Cephalaspidomorphi"
    if "Myxini" in p: return "Myxini"
    return ncbi_class or None


def choose(obs: list[Obs], group: str) -> Obs | None:
    """The chosen maximum lifespan of one species from its matched observations (their order breaks ties)."""
    with_value = [o for o in obs if o.value_years is not None]
    anage = next((o for o in with_value if o.source == "AnAge"), None)
    if anage: return anage
    for src in PRECEDENCE[group]:
        if src == "AnAge": continue
        h = [o for o in with_value if o.source == src]
        if src == "SeaLifeBase": h = [o for o in h if o.aquatic is True]
        if not h: continue
        for sel in (lambda o: o.record_type != "estimated" and o.quality != "questionable",
                    lambda o: o.record_type != "estimated",
                    lambda o: o.record_type == "estimated"):
            hh = [o for o in h if sel(o)]
            if not hh: continue
            best = hh[0]
            for o in hh:
                tie_with_origin = o.value_years == best.value_years and best.record_type == "observed" and o.record_type != "observed"
                if o.value_years > best.value_years or tie_with_origin: best = o
            return best
    return None


def load_all(b2_in: dict) -> tuple[list[AnAgeEntry], list[Obs], dict]:
    """AnAge entries, every observation in the order the rules read them, and the source synonym tables."""
    anage = load_anage(b2_in["anage"])
    fb, fb_syn = load_fishbase_like(b2_in["fishbase"], "fb")
    slb, slb_syn = load_fishbase_like(b2_in["sealifebase"], "slb")
    usgs = load_usgs(b2_in["usgs"])
    from . import euring
    eur = euring.load(b2_in["euring"])
    obs = anage_observations(anage) + fb + slb + eur + usgs
    return anage, obs, {"FishBase": fb_syn, "SeaLifeBase": slb_syn}


def matched(o: Obs, m: Match) -> Obs:
    return replace(o, species_taxid=m.species_taxid, match_method=m.method)
