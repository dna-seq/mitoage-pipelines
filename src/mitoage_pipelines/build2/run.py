"""The stages of `mitoage build2`, each writing its table, and the choice of inputs (pinned or fresh)."""
import json
import time
from pathlib import Path

import polars as pl

from .. import paths
from . import assemble as asm
from . import pins

EFETCH_LABEL = "build2-653afe2a3b"


def label(folder: Path, default: str) -> str:
    p = folder / "label.json"
    return json.loads(p.read_text())["label"] if p.is_file() else default


def release_file(release: str) -> tuple[Path, str]:
    """`237` (or another number): the release file; `efetch`: build 2's versions fetched by E-utilities."""
    if release == "efetch":
        p = paths.B2_REFSEQ / f"efetch-{EFETCH_LABEL}" / "mitochondrion.efetch.gbff.gz"
        return p, label(p.parent, "build 2's RefSeq accession.versions by E-utilities efetch")
    p = paths.refseq_release(release)
    return p, label(p.parent, pins.REFSEQ_VERSION if release == pins.REFSEQ_RELEASE else f"release {release}")


def taxdump_folder(name: str) -> tuple[Path, str]:
    p = paths.taxdump(name)
    if name == f"current-{pins.TAXDUMP_BUILD2}": return p, pins.TAXONOMY_VERSION
    return p, label(p, f"new_taxdump archive of {name}" if not name.startswith("current-") else f"new_taxdump of {name[8:]}")


def inputs(release: str, taxonomy: str) -> asm.Inputs:
    rf, rl = release_file(release)
    tf, tl = taxdump_folder(taxonomy)
    for p in (rf, tf):
        if not p.exists():
            raise FileNotFoundError(f"no {p}: download it first (`mitoage download refseq` / `mitoage download taxonomy`)")
    return asm.Inputs(build1=paths.B2_BUILD1, release=rf, release_label=rl, taxdump=tf, taxonomy_label=tl, anage=paths.B2_ANAGE,
                      fishbase=paths.B2_FISHBASE / pins.FISHBASE_VERSION, sealifebase=paths.B2_SEALIFEBASE / pins.FISHBASE_VERSION,
                      usgs=paths.B2_USGS / "longevity.csv", euring=paths.B2_EURING, genbank=paths.B2_GENBANK, legacy=paths.B2_LEGACY)


def state(inp: asm.Inputs, sources: bool = True) -> asm.State:
    st = asm.State(inp)
    asm.load_build1(st)
    if sources: asm.load_sources(st)
    return st


def stage_refseq(inp, interim: Path) -> None:
    st = state(inp)
    choice = asm.one_per_species(st.records)
    chosen = {recs[0].header.version for recs in choice.values()}
    df = pl.DataFrame([{"version": r.header.version, "length": r.header.length, "update_date": r.header.update_date,
                        "definition": r.header.definition, "taxid": r.header.taxid, "taxid_current": r.taxid_current,
                        "species_taxid": r.species_taxid, "species_class": r.species_class, "completeness": r.completeness,
                        "length_qc": r.length_qc, "has_sequence": r.header.has_sequence, "chosen": r.header.version in chosen}
                       for r in st.records])
    interim.mkdir(parents=True, exist_ok=True)
    df.write_parquet(interim / "refseq_records.parquet")
    print(f"{df.height} metazoan records, {df['completeness'].str.starts_with('complete').sum()} complete, "
          f"{len(choice)} species with a complete mitogenome -> {interim / 'refseq_records.parquet'}")


def stage_match(inp, interim: Path) -> None:
    st = state(inp)
    rows = []
    for o in st.obs:
        m = st.matches[asm.ls.rec_key(o.source, o.source_name, o.speccode)]
        rows.append({"source": o.source, "source_record_id": o.source_record_id, "source_name": o.source_name, "value_years": o.value_years,
                     "record_type": o.record_type, "quality": o.quality, "match_method": m.method, "match_status": m.status,
                     "species_taxid": m.species_taxid, "candidates": ";".join(map(str, m.candidates))})
    df = pl.DataFrame(rows, infer_schema_length=None)
    interim.mkdir(parents=True, exist_ok=True)
    df.write_parquet(interim / "observations_matched.parquet")
    print(df.group_by("source", "match_status").len().sort("source", "match_status"))


def stage_status(inp, interim: Path) -> dict:
    st = state(inp)
    pl_ = asm.plan(st)
    rows = [{"hagrid": o["id"], "scientific_name": o["s"].species, "cited": o["s"].refseq, "version": o["version"],
             "change": o["change"][2] if o["change"] else None, "retired_reason": o["retired"], "ncbi_species_taxid": o["eff"],
             "lifespan_source": o["lifespan"].source if o["lifespan"] else None,
             "max_longevity_years": o["lifespan"].value_years if o["lifespan"] else o["s"].max_longevity} for o in pl_["origs"]]
    cands = [{"ncbi_species_taxid": c["t"], "refseq": c["rec"].header.version, "group": c["group"], "lifespan_source": c["chosen"].source,
              "max_longevity_years": c["chosen"].value_years} for c in pl_["cands"]]
    interim.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_parquet(interim / "build1_status.parquet")
    pl.DataFrame(cands).write_parquet(interim / "candidates.parquet")
    df = pl.DataFrame(rows)
    print(f"build 1 species: {df.height}; same version {(df['change'].is_null() & df['retired_reason'].is_null()).sum()}, "
          f"new version {(df['change'] == 'upstream_version_change').sum()}, successor {(df['change'] == 'refseq_successor').sum()}, "
          f"retired {df['retired_reason'].is_not_null().sum()}; candidates (new species): {len(cands)}")
    return {"state": st, "plan": pl_}


def stage_assemble(inp, out: Path) -> dict:
    from . import output
    t0 = time.time()
    st = state(inp)
    pl_ = asm.plan(st)
    b2 = asm.assemble(st, pl_)
    info = {"build": "2 (reproduction)", "build_date": inp.build_date, "refseq": inp.release_label, "refseq_file": str(inp.release),
            "taxonomy": inp.taxonomy_label, "taxonomy_folder": str(inp.taxdump), "anage_build": st.anage_build,
            "seconds": round(time.time() - t0, 1)}
    n = output.write_all(b2, out, info)
    status = {}
    for s in b2["species"]: status[s["record_status"]] = status.get(s["record_status"], 0) + 1
    print(f"{n} -> {out}\nrecord status: {status}; {b2['report']}")
    return b2


def stage_primary_list(inp, interim: Path) -> Path:
    """The build 1 species whose INSDC source record the re-annotation test needs: those that keep their
    accession.version and whose recomputed region or codon values differ from build 1 (assembled once
    without any INSDC record). Build 2's audit chose its 160 records by hand; this is the rule written down."""
    st = state(inp)
    pl_ = asm.plan(st)
    b2 = asm.assemble(st, pl_, use_primary=False)
    same = {o["id"] for o in pl_["origs"] if not o["change"] and not o["retired"]}
    want = sorted({int(c["row_key"].split(":")[0]) for c in b2["changelog"]
                   if c["table"] in ("region_stats", "codon_usage") and c["field"] != "(row)" and c["row_key"] != "*"} & same)
    interim.mkdir(parents=True, exist_ok=True)
    out = interim / "primary_hagrids.txt"
    out.write_text("".join(f"{h}\n" for h in want))
    print(f"{len(want)} build 1 species need their INSDC source record -> {out}")
    return out
