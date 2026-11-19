"""H27: the species sets an analysis can be restricted to, by the confidence level of a record and by how well
its lifespan is sampled.

    not_unreliable   the confidence level is not "unreliable"
    high             the level is "high"
    well_sampled     not_unreliable, and the lifespan is AnAge's with a sample size of medium, large or huge
    captive, wild    not_unreliable, and the lifespan record is of that kind

A set is chosen for a whole run by the environment variable MITOAGE_SPECIES_SET, so that an analysis changes in
nothing but its species: `aa_pgls.default_dataset`, which every analysis of lifespan calls, then keeps only the
species of the set. Unset, nothing is restricted and the analyses are those of their own entries.

The level is read from `data_confidence` in the build's species table. A table without that column (a build made
before the level existed) is refused: a missing level is never taken as "high" or as "not unreliable".
"""
import os
from pathlib import Path

import pandas as pd

ENV = "MITOAGE_SPECIES_SET"
SETS = ("not_unreliable", "high", "well_sampled", "captive", "wild")
LEVELS = ("high", "caveat", "low", "unreliable")
WELL_SAMPLED = ("medium", "large", "huge")
COLUMNS = ["slug", "data_confidence", "lifespan_source", "lifespan_sample_size", "lifespan_record_type"]


def active() -> str | None:
    """The set named by MITOAGE_SPECIES_SET, None when it is unset or empty; an unknown name is refused."""
    name = os.environ.get(ENV, "").strip()
    if not name:
        return None
    if name not in SETS:
        raise SystemExit(f"{ENV}={name!r}: not a species set; the sets are {', '.join(SETS)}")
    return name


def members(species: pd.DataFrame, name: str) -> pd.Series:
    """True for the rows of `species` (the build's species table) in the set `name`."""
    if name not in SETS:
        raise ValueError(f"{name!r} is not a species set; the sets are {', '.join(SETS)}")
    if "data_confidence" not in species.columns:
        raise SystemExit("the species table has no data_confidence column: this build gives no confidence level, "
                         "so no restricted species set can be made from it (download a build that has one)")
    level = species.data_confidence
    unknown = sorted(set(level.dropna()) - set(LEVELS))
    if unknown:
        raise SystemExit(f"data_confidence holds values that are no level: {', '.join(map(str, unknown))}")
    ok = level.isin(LEVELS) & level.ne("unreliable")           # a missing level is in no set
    if name == "not_unreliable":
        return ok
    if name == "high":
        return level.eq("high")
    if name == "well_sampled":
        return ok & species.lifespan_source.eq("AnAge") & species.lifespan_sample_size.isin(WELL_SAMPLED)
    return ok & species.lifespan_record_type.eq(name)           # captive, wild


def slugs(build: Path, name: str) -> set[str]:
    """The slugs of the set `name`, from the species table of the build in `build`."""
    species = pd.read_parquet(build / "species.parquet")
    missing = [c for c in COLUMNS if c not in species.columns and c != "data_confidence"]
    if missing:
        raise SystemExit(f"the species table lacks {', '.join(missing)}")
    return set(species.slug[members(species, name)])
