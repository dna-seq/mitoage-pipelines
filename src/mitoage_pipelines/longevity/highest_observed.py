"""H27: a second lifespan per species, the highest observed value across the build's lifespan sources.

The build takes AnAge's value whenever AnAge has one. A maximum lifespan is a record, so a documented higher
value from another source is the better estimate unless it is an error. This writes a copy of a build in which
each species' lifespan is the largest of its qualifying observations, so that any analysis run on the copy is
changed in nothing but the lifespan.

A qualifying observation: an observed age (record type captive, wild or observed; not estimated, not the 2015
database's compiled value), rated high or acceptable by its source's scale. The species' chosen value always
stays in the running, whatever its kind, so no species loses its lifespan and none gets a lower one.
"""
from pathlib import Path

import pandas as pd

OBSERVED = ("captive", "wild", "observed")
QUALITY_OK = ("high", "acceptable")
COLUMNS = {"value_years": "max_longevity_years", "source": "lifespan_source", "source_record_id": "lifespan_source_record_id",
           "record_type": "lifespan_record_type", "quality": "lifespan_quality"}


def highest(species: pd.DataFrame, observations: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """-> (the species table with the highest observed lifespan, the species whose value changed: slug, class,
    the two values, their sources and the ratio). Ties keep the chosen value."""
    o = observations
    ok = o.record_type.isin(OBSERVED) & o.quality.isin(QUALITY_OK)     # the chosen value may qualify: it only ties itself
    best = (o[ok].sort_values(["species_id", "value_years", "source", "source_record_id"], ascending=[True, False, True, True])
            .drop_duplicates("species_id").set_index("species_id"))
    out = species.copy()
    cand = out.id.map(best.value_years)
    higher = cand.gt(out.max_longevity_years)                   # NaN (no candidate) is False; a tie keeps the chosen value
    changed = pd.DataFrame({
        "slug": out.slug[higher], "class": out["class"][higher],
        "chosen_years": out.max_longevity_years[higher], "chosen_source": out.lifespan_source[higher],
        "highest_years": cand[higher], "highest_source": out.id[higher].map(best.source),
    })
    changed["ratio"] = changed.highest_years / changed.chosen_years
    for src, dst in COLUMNS.items():
        out.loc[higher, dst] = out.id[higher].map(best[src])
    if "lifespan_sample_size" in out.columns:                    # AnAge's sample size is not the new record's
        out.loc[higher & out.lifespan_source.ne("AnAge"), "lifespan_sample_size"] = None
    return out, changed.sort_values("ratio", ascending=False).reset_index(drop=True)


def main(build: Path, dest: Path, table: Path) -> pd.DataFrame:
    """Writes `dest`, a build folder like `build` whose species table has the highest observed lifespan (every
    other file is a link to the original), and `table`, the species whose value changed."""
    if dest.resolve() == build.resolve():
        raise SystemExit("the copy must not be written over the build it is made from")
    species = pd.read_parquet(build / "species.parquet")
    out, changed = highest(species, pd.read_parquet(build / "lifespan_observations.parquet"))
    dest.mkdir(parents=True, exist_ok=True)
    for f in sorted(build.iterdir()):
        link = dest / f.name
        if f.name == "species.parquet" or link.exists() or link.is_symlink():
            continue
        link.symlink_to(f.resolve(), target_is_directory=f.is_dir())
    out.to_parquet(dest / "species.parquet", index=False)
    table.parent.mkdir(parents=True, exist_ok=True)
    changed.to_csv(table, sep="\t", index=False, float_format="%.6g")
    print(f"{len(changed)} of {len(species)} species get a higher lifespan; by class and source:")
    print(changed.groupby(["class", "highest_source"]).size().to_string())
    print(changed.ratio.describe().to_string())
    return changed
