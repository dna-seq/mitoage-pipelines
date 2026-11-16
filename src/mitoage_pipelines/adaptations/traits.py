"""The trait states of H21-H25 for every build species in each trait's scope (docs/HYPOTHESES.md).

One row per (trait, species in scope): state 1, 0, or unknown (NA), and how it was decided. A species outside a
trait's scope has no row. Unknown is never written as 0.

- diadromy (H21; ray-finned fishes): FishBase v26.06 `species.AnaCat` anadromous, catadromous, amphidromous or
  diadromous = 1; potamodromous, oceanodromous, non-migratory or oceano-estuarine = 0; empty, or no FishBase match,
  unknown. migration (H21, secondary): any of the six migratory values = 1, non-migratory = 0, oceano-estuarine,
  empty or no match unknown.
- deep_habitat (H22; ray-finned fishes): `species.DemersPelag` bathydemersal or bathypelagic = 1, any other value = 0,
  empty or no match unknown. deep_depth (H22, sensitivity): `species.DepthRangeDeep` >= 1000 m = 1, < 200 m = 0,
  200 to 999 m and empty unknown ("intermediate" is said as such).
- marine (H23; mammals): cetacean families (under Cetacea or Artiodactyla), Otariidae, Phocidae, Odobenidae, Sirenia,
  Enhydra = 1, other mammals 0.
- hibernation (H23, H24; mammals): Heldstab et al. 2018 with Wilkinson & Adams 2019 for bats (`hibernation.combine`).
- bat_hibernation (H25; Chiroptera): Wilkinson & Adams 2019 alone (`hibernation.match`).
- bats (H23's positive control; mammals): order Chiroptera = 1, other mammals 0.

FishBase names are matched by H20's rule (`thermal.fishbase_codes`).
"""
from pathlib import Path

import pandas as pd

from ..longevity.thermal import fishbase_codes
from . import hibernation as hib

DIADROMOUS = ("anadromous", "catadromous", "amphidromous", "diadromous")
OTHER_ANACAT = ("potamodromous", "oceanodromous", "non-migratory", "oceano-estuarine")
MIGRATORY = DIADROMOUS + ("potamodromous", "oceanodromous")
BATHY = ("bathydemersal", "bathypelagic")
DEEP_M, SHALLOW_M = 1000.0, 200.0
CETACEAN_FAMILIES = ("Balaenidae", "Balaenopteridae", "Delphinidae", "Eschrichtiidae", "Hyperoodontidae", "Iniidae",
                     "Kogiidae", "Lipotidae", "Monodontidae", "Phocoenidae", "Physeteridae", "Ziphiidae", "Platanistidae")
MARINE_FAMILIES = CETACEAN_FAMILIES + ("Otariidae", "Phocidae", "Odobenidae")
SCOPE = {"diadromy": "Actinopterygii", "migration": "Actinopterygii", "deep_habitat": "Actinopterygii",
         "deep_depth": "Actinopterygii", "marine": "Mammalia", "hibernation": "Mammalia",
         "bat_hibernation": "Mammalia", "bats": "Mammalia"}
COLUMNS = ["trait", "slug", "class", "order", "family", "state", "how"]


def _rows(trait: str, sp: pd.DataFrame, state: pd.Series, how: pd.Series) -> pd.DataFrame:
    return pd.DataFrame({"trait": trait, "slug": sp.slug.to_numpy(), "class": sp["class"].to_numpy(),
                         "order": sp.order.to_numpy(), "family": sp.family.to_numpy(),
                         "state": pd.array(state.to_numpy(), dtype="Int64"), "how": how.to_numpy()})


def fish_traits(fish: pd.DataFrame, fb: Path) -> list[pd.DataFrame]:
    codes = fishbase_codes(fb, fish).set_index("slug")
    fbs = pd.read_parquet(fb / "species.parquet", columns=["SpecCode", "AnaCat", "DemersPelag", "DepthRangeDeep"]
                          ).drop_duplicates("SpecCode").set_index("SpecCode")
    code = fish.slug.map(codes.SpecCode)
    matched = code.notna()
    get = lambda col: code.map(fbs[col]).where(matched)
    ana, dem, depth = get("AnaCat"), get("DemersPelag"), get("DepthRangeDeep").astype(float)
    nomatch = "no FishBase match: " + fish.slug.map(codes.how).fillna("")
    out = []
    s = pd.Series(pd.NA, index=fish.index, dtype="Int64")
    st = s.copy(); st[ana.isin(DIADROMOUS)] = 1; st[ana.isin(OTHER_ANACAT)] = 0
    how = pd.Series("AnaCat " + ana.astype(str), index=fish.index).where(ana.notna(), "AnaCat empty").where(matched, nomatch)
    out.append(_rows("diadromy", fish, st, how))
    st = s.copy(); st[ana.isin(MIGRATORY)] = 1; st[ana.eq("non-migratory")] = 0
    out.append(_rows("migration", fish, st, how))
    st = s.copy(); st[dem.isin(BATHY)] = 1; st[dem.notna() & ~dem.isin(BATHY)] = 0
    how = pd.Series("DemersPelag " + dem.astype(str), index=fish.index).where(dem.notna(), "DemersPelag empty").where(matched, nomatch)
    out.append(_rows("deep_habitat", fish, st, how))
    st = s.copy(); st[depth >= DEEP_M] = 1; st[depth < SHALLOW_M] = 0
    how = pd.Series("DepthRangeDeep " + depth.astype(str), index=fish.index)
    how[depth.isna()] = "DepthRangeDeep empty"
    how[depth.ge(SHALLOW_M) & depth.lt(DEEP_M)] = "DepthRangeDeep 200 to 999 m: intermediate, left out"
    out.append(_rows("deep_depth", fish, st, how.where(matched, nomatch)))
    return out


def mammal_traits(mam: pd.DataFrame, hib_dir: Path, names_dmp: Path | None) -> list[pd.DataFrame]:
    out = []
    marine = mam.family.isin(MARINE_FAMILIES) | mam.order.isin(["Cetacea", "Sirenia"]) | mam.genus.eq("Enhydra")
    out.append(_rows("marine", mam, marine.astype(int), pd.Series("taxonomic rule", index=mam.index)))
    bat = mam.order.eq("Chiroptera")
    out.append(_rows("bats", mam, bat.astype(int), pd.Series("order", index=mam.index)))
    bats_file, mammals_file = hib_dir / hib.BATS_FILE, hib_dir / hib.MAMMALS_FILE
    bats_t = hib.read_bats(bats_file) if bats_file.is_file() else None
    mam_t = hib.read_mammals(mammals_file) if mammals_file.is_file() else None
    if bats_t is None and mam_t is None:
        print(f"NOTE: no hibernation table under {hib_dir}: hibernation is unknown for every mammal (`download hibernation`)")
    c = hib.combine(mam, mam_t, bats_t, names_dmp).set_index("slug")
    out.append(_rows("hibernation", mam, mam.slug.map(c.state), mam.slug.map(c.how)))
    chiro = mam[bat]
    if bats_t is not None:
        m = hib.match(chiro, bats_t, names_dmp).set_index("slug")
        out.append(_rows("bat_hibernation", chiro, chiro.slug.map(m.state), chiro.slug.map(m.how)))
    else:
        out.append(_rows("bat_hibernation", chiro, pd.Series(pd.NA, index=chiro.index),
                         pd.Series("no bat table", index=chiro.index)))
    return out


def classify(species: pd.DataFrame, fb: Path, hib_dir: Path, names_dmp: Path | None) -> pd.DataFrame:
    sp = species[species.record_status.ne("retired")].reset_index(drop=True)
    parts = fish_traits(sp[sp["class"].eq("Actinopterygii")].reset_index(drop=True), fb)
    parts += mammal_traits(sp[sp["class"].eq("Mammalia")].reset_index(drop=True), hib_dir, names_dmp)
    t = pd.concat(parts, ignore_index=True)[COLUMNS]
    return t.sort_values(["trait", "slug"], kind="stable").reset_index(drop=True)


def main(data: Path, fb: Path, hib_dir: Path, names_dmp: Path | None, out: Path) -> pd.DataFrame:
    t = classify(pd.read_parquet(data / "species.parquet"), fb, hib_dir, names_dmp)
    out.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(out, sep="\t", index=False)
    summary = t.assign(state=t.state.astype("string").fillna("unknown")).groupby(["trait", "state"]).size().unstack(fill_value=0)
    print(summary.to_string())
    return t
