"""The build 2 tables as parquet, with the columns and types of the published ones (site repository:
pipeline/steps/output.ts), plus review_list.csv and build-info.json."""
import csv
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .genome import AMINO_ACIDS, CODONS

PROV = [("added_in_build", pa.int32()), ("added_at", pa.string()), ("last_modified_build", pa.int32()),
        ("last_modified_at", pa.string()), ("record_status", pa.string()), ("source_version", pa.string())]
TRAIT_COLS = ["adult_weight_g", "birth_weight_g", "metabolic_rate_w", "body_temperature_k", "female_maturity_days",
              "male_maturity_days", "gestation_days", "litter_size", "litters_per_year", "growth_rate_per_day"]
I, S, D, B = pa.int32(), pa.string(), pa.float64(), pa.bool_()


def write(path: Path, cols: list[tuple[str, pa.DataType]], rows: list[dict]) -> int:
    table = pa.table({n: pa.array([r.get(n) for r in rows], type=t) for n, t in cols})
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path)
    return len(rows)


def species_rows(species: list[dict], extra: dict) -> list[dict]:
    out = []
    for s in species:
        e, l, t = extra[s["id"]], s["lifespan"], s["traits"]
        out.append({**{k: s[k] for k in ("id", "slug", "scientific_name", "common_name", "refseq", "kingdom", "phylum", "class", "order",
                                          "family", "genus", "max_longevity_years", "mtdna_size_bp", "non_coding_bp", "gc_percent", "legacy_id",
                                          "ncbi_taxid", "ncbi_species_taxid", "retired_reason", "added_in_build", "added_at",
                                          "last_modified_build", "last_modified_at", "record_status", "source_version")},
                    **{c: t.get(c) for c in TRAIT_COLS}, "traits_source": t.get("source"),
                    "anage_build15_hagrid": e["anage_hagrid"], "anage_match": e["anage_match"],
                    "anage_build15_max_longevity_years": e["anage_max"], "transl_table": e["transl_table"],
                    "lifespan_source": l["source"], "lifespan_source_version": l["source_version"],
                    "lifespan_source_record_id": l["source_record_id"], "lifespan_quality": l["quality"],
                    "lifespan_sample_size": l["sample_size"], "lifespan_record_type": l["record_type"], "lifespan_license": l["license"],
                    "previous_scientific_names": "; ".join(s["previous"]) or None, "review_flags": "; ".join(s["review_flags"]) or None,
                    "change_kinds": "; ".join(s["change_kinds"]) or None})
    return out


SPECIES_COLS = ([("id", I), ("slug", S), ("scientific_name", S), ("common_name", S), ("refseq", S), ("kingdom", S), ("phylum", S),
                 ("class", S), ("order", S), ("family", S), ("genus", S), ("max_longevity_years", D), ("mtdna_size_bp", I),
                 ("non_coding_bp", I), ("gc_percent", D), ("legacy_id", I)] + [(c, D) for c in TRAIT_COLS] +
                [("traits_source", S), ("anage_build15_hagrid", S), ("anage_match", S), ("anage_build15_max_longevity_years", D),
                 ("ncbi_taxid", I), ("ncbi_species_taxid", I), ("transl_table", I), ("lifespan_source", S), ("lifespan_source_version", S),
                 ("lifespan_source_record_id", S), ("lifespan_quality", S), ("lifespan_sample_size", S), ("lifespan_record_type", S),
                 ("lifespan_license", S), ("previous_scientific_names", S), ("retired_reason", S), ("review_flags", S), ("change_kinds", S)] + PROV)
REGION_COLS = [("species_id", I), ("slug", S), ("region", S), ("length_bp", I), ("a", I), ("c", I), ("g", I), ("t", I), ("other", I),
               ("gc_percent", D), ("gc_skew", D), ("at_skew", D), ("source", S)] + PROV
CODON_COLS = ([("species_id", I), ("slug", S), ("region", S), ("transl_table", I), ("source", S)] + [(c, I) for c in CODONS] +
              [(f"aa_{a}", I) for a in AMINO_ACIDS] + [("aa_stop", I)] + PROV)
FEATURE_COLS = [("species_id", I), ("slug", S), ("feature_key", S), ("name", S), ("type", S), ("start", I), ("end", I), ("strand", I)] + PROV
OBS_COLS = [("species_id", I), ("slug", S), ("ncbi_species_taxid", I), ("source", S), ("source_record_id", S), ("source_name", S),
            ("value_years", D), ("record_type", S), ("quality", S), ("sample_size_text", S), ("sex", S), ("detail", S), ("reference", S),
            ("license", S), ("match_method", S), ("is_chosen", B)] + PROV
CHANGELOG_COLS = [("build", I), ("changed_at", S), ("table", S), ("row_key", S), ("field", S), ("old_value", S), ("new_value", S),
                  ("change_type", S), ("diff_class", S), ("reason", S), ("upstream_source", S), ("upstream_version", S), ("scope", S)]
REVIEW_COLS = ["review_type", "ncbi_taxid", "scientific_name", "refseq", "source", "source_name", "detail"]


def write_all(b2: dict, out: Path, info: dict) -> dict:
    slug = {s["id"]: s["slug"] for s in b2["species"]}
    with_slug = lambda rows: [{**r, "slug": slug.get(r["species_id"])} for r in rows]
    usage = [{**{k: v for k, v in u.items() if k not in ("counts", "amino_acids")}, **u["counts"],
              **{f"aa_{a}": u["amino_acids"][a] for a in AMINO_ACIDS}, "aa_stop": u["amino_acids"]["*"]} for u in b2["usage"]]
    n = {"species": write(out / "species.parquet", SPECIES_COLS, species_rows(b2["species"], b2["extra"])),
         "region_stats": write(out / "region_stats.parquet", REGION_COLS, with_slug(b2["stats"])),
         "codon_usage": write(out / "codon_usage.parquet", CODON_COLS, with_slug(usage)),
         "genome_features": write(out / "genome_features.parquet", FEATURE_COLS, with_slug(b2["features"])),
         "lifespan_observations": write(out / "lifespan_observations.parquet", OBS_COLS, with_slug(b2["observations"])),
         "changelog": write(out / "changelog.parquet", CHANGELOG_COLS, b2["changelog"])}
    with open(out / "review_list.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n"); w.writerow(REVIEW_COLS)
        w.writerows([[r[c] for c in REVIEW_COLS] for r in b2["review"]])
    n["review_list"] = len(b2["review"])
    (out / "build-info.json").write_text(json.dumps({**info, "rows": n, "report": b2["report"]}, indent=1, default=str) + "\n")
    return n
