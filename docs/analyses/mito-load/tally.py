"""Design tally for docs/MITO_LOAD_ADAPTATIONS.md: how many species of the build carry each candidate adaptation,
how many are unknown, how often the trait changes on the dated tree, and how many trait-bearing species have
nuclear partners. It is design information, read before H21-H23 were written; it tests nothing.

    uv run python docs/analyses/mito-load/tally.py SITE_CHECKOUT FISHBASE_DIR OUT_TSV [BATS_TABLE NAMES_DMP [MAMMAL_TABLE]]

SITE_CHECKOUT: dna-seq/mitoage at 03c0949 or later with the same data (public/data/species.csv,
src/data/generated/phylogeny.json, public/data/mitoport/partners_species.parquet).
FISHBASE_DIR: FishBase v26.06 from source.coop, the species, synonyms and stocks tables (`download fishbase-temperature`).
BATS_TABLE, NAMES_DMP: the bat table of `download hibernation` and NCBI's names.dmp; with them, bat hibernation is tallied.
MAMMAL_TABLE: Heldstab et al. 2018's sheet; with it, hibernation across mammals (both sources combined) is tallied.

Every rule below is an explicit taxonomic rule or a FishBase field. A trait whose classification needs a source we
do not hold (hibernation, high altitude, migration in birds, diving in birds) is not tallied here.
"""
import json
import sys
from pathlib import Path

import pandas as pd

from mitoage_pipelines.adaptations import hibernation as hib
from mitoage_pipelines.longevity.thermal import fishbase_codes

FISHES = ("Actinopterygii", "Chondrichthyes")
CETACEAN_FAMILIES = {"Balaenidae", "Balaenopteridae", "Delphinidae", "Eschrichtiidae", "Hyperoodontidae", "Iniidae",
                     "Kogiidae", "Lipotidae", "Monodontidae", "Phocoenidae", "Physeteridae", "Ziphiidae", "Platanistidae"}


def isin(col, values):
    return lambda d: d[col].isin(values)


def rules():
    """name -> (classes in scope, is_trait(d) -> bool Series, is_unknown(d) -> bool Series or None)"""
    no = None
    air = {"Facultative", "Obligate", "FacultativeGenus", "FacultativeObligate"}
    return {
        "mammal flight: order Chiroptera": (("Mammalia",), isin("order", {"Chiroptera"}), no),
        "mammal marine: cetacean families (under Cetacea or Artiodactyla), Otariidae, Phocidae, Odobenidae, Sirenia, Enhydra":
            (("Mammalia",), lambda d: d.family.isin(CETACEAN_FAMILIES | {"Otariidae", "Phocidae", "Odobenidae"})
             | d.order.isin({"Cetacea", "Sirenia"}) | d.genus.eq("Enhydra"), no),
        "mammal semi-aquatic: Lutrinae genera, Castor, Ornithorhynchus, Neomys, Myocastor, Hippopotamidae, Condylura":
            (("Mammalia",), lambda d: d.genus.isin({"Aonyx", "Hydrictis", "Lontra", "Lutra", "Pteronura", "Castor",
                                                   "Ornithorhynchus", "Neomys", "Myocastor", "Condylura"})
             | d.family.eq("Hippopotamidae"), no),
        "mammal subterranean: Bathyergidae, Spalacidae, Talpidae, Spalacopus, Ellobius":
            (("Mammalia",), lambda d: d.family.isin({"Bathyergidae", "Spalacidae", "Talpidae"})
             | d.genus.isin({"Spalacopus", "Ellobius"}), no),
        "mammal myrmecophagy: Myrmecophagidae, Tachyglossidae, Myrmecobiidae, Orycteropodidae, Proteles, Priodontes, Tolypeutes, Melursus":
            (("Mammalia",), lambda d: d.family.isin({"Myrmecophagidae", "Tachyglossidae", "Myrmecobiidae", "Orycteropodidae"})
             | d.genus.isin({"Proteles", "Priodontes", "Tolypeutes", "Melursus"}), no),
        "bird flightless: Struthioniformes, Casuariiformes, Sphenisciformes, Strigops":
            (("Aves",), lambda d: d.order.isin({"Struthioniformes", "Casuariiformes", "Sphenisciformes"})
             | d.genus.eq("Strigops"), no),
        "fish diadromous: FishBase AnaCat anadromous, catadromous, amphidromous or diadromous":
            (("Actinopterygii",), isin("AnaCat", {"anadromous", "catadromous", "amphidromous", "diadromous"}),
             lambda d: d.AnaCat.isna()),
        "fish deep: FishBase DepthRangeDeep >= 1000 m":
            (("Actinopterygii",), lambda d: d.DepthRangeDeep >= 1000, lambda d: d.DepthRangeDeep.isna()),
        "fish deep: FishBase DemersPelag bathydemersal or bathypelagic":
            (("Actinopterygii",), isin("DemersPelag", {"bathydemersal", "bathypelagic"}), lambda d: d.DemersPelag.isna()),
        "fish polar: FishBase EnvTemp 'polar' in any stock":
            (("Actinopterygii",), lambda d: d.EnvTemp.fillna("").str.contains("polar"), lambda d: d.EnvTemp.isna()),
        "fish air-breathing: FishBase AirBreathing facultative or obligate ('WaterAssumed' = unknown)":
            (("Actinopterygii",), isin("AirBreathing", air), lambda d: d.AirBreathing.isna() | d.AirBreathing.eq("WaterAssumed")),
        "fish regional endothermy: Thunnus, Katsuwonus, Istiophoridae, Xiphiidae, Lampris":
            (("Actinopterygii",), lambda d: d.genus.isin({"Thunnus", "Katsuwonus", "Lampris"})
             | d.family.isin({"Istiophoridae", "Xiphiidae"}), no),
        "shark regional endothermy: Lamnidae, Alopiidae":
            (("Chondrichthyes",), isin("family", {"Lamnidae", "Alopiidae"}), no),
        "fish electric discharge: FishBase Electrogenic weakly or strongly discharging":
            (FISHES, lambda d: d.Electrogenic.fillna("").str.contains("discharging"), lambda d: d.Electrogenic.isna()),
    }


def min_changes(tree: dict, state: dict) -> int | None:
    """Fewest changes of a 0/1 trait on the tree (Sankoff, unit costs; polytomies allowed). Tips in `state` with
    value None are unknown and cost nothing either way; tips not in `state` are pruned."""
    INF = 10 ** 9
    stack, post = [tree], []
    while stack:                       # iterative post-order: the tree is deep
        n = stack.pop(); post.append(n); stack.extend(n.get("children", []))
    cost = {}
    for n in reversed(post):
        if "leaf" in n:
            if n["leaf"] not in state: cost[id(n)] = None; continue
            v = state[n["leaf"]]
            cost[id(n)] = (0, 0) if v is None else ((0, INF) if v == 0 else (INF, 0))
            continue
        cs = [cost[id(c)] for c in n["children"] if cost[id(c)] is not None]
        cost[id(n)] = None if not cs else tuple(sum(min(c[0] + (k != 0), c[1] + (k != 1)) for c in cs) for k in (0, 1))
    r = cost[id(tree)]
    return None if r is None else min(r)


def trait_clades(tree: dict, state: dict) -> int:
    """Maximal clades whose known tips (at least one) are all trait-bearing. Each needs a gain on its stem unless
    an ancestor gained the trait and its other descendants lost it, so it is not a bound on gains; with
    `min_changes` it brackets how many times the trait arose within the build's sample."""
    out, stack, post = 0, [tree], []
    while stack:
        n = stack.pop(); post.append(n); stack.extend(n.get("children", []))
    pure = {}                          # id -> None (no known tip), True (all known are 1), False
    for n in reversed(post):
        if "leaf" in n:
            v = state.get(n["leaf"]); pure[id(n)] = None if v is None else (v == 1); continue
        vals = [pure[id(c)] for c in n["children"] if pure[id(c)] is not None]
        pure[id(n)] = None if not vals else all(vals)
    def count(n, parent_pure):
        nonlocal out
        p = pure[id(n)]
        if p and not parent_pure: out += 1
        for c in n.get("children", []): count(c, bool(p))
    sys.setrecursionlimit(100_000); count(tree, False)
    return out


def main(site: Path, fb: Path, out: Path, bats: Path | None = None, names: Path | None = None,
         mammals: Path | None = None) -> pd.DataFrame:
    sp = pd.read_csv(site / "public/data/species.csv", low_memory=False)
    sp = sp[sp.retired_reason.isna()].copy()
    fish = sp["class"].isin(FISHES)
    match = fishbase_codes(fb, sp[fish])
    fbs = pd.read_parquet(fb / "species.parquet", columns=["SpecCode", "AnaCat", "DepthRangeDeep", "DemersPelag",
                                                           "AirBreathing", "Electrogenic"])
    st = pd.read_parquet(fb / "stocks.parquet", columns=["SpecCode", "EnvTemp"]).dropna()
    env = st.groupby("SpecCode").EnvTemp.agg(lambda s: "|".join(sorted(set(s)))).rename("EnvTemp")
    match = match.merge(fbs, on="SpecCode", how="left").merge(env, left_on="SpecCode", right_index=True, how="left")
    sp = sp.merge(match.drop(columns="SpecCode"), on="slug", how="left")
    r = rules()
    if bats is not None:
        chiro = sp[sp.order.eq("Chiroptera")]
        m = hib.match(chiro, hib.read_bats(bats), names).set_index("slug")
        sp["bat_hibernation"] = sp.slug.map(m.state)
        r["bat hibernation: Wilkinson & Adams 2019 (DRUM), within Chiroptera"] = (
            ("Mammalia",), lambda d: d.bat_hibernation.eq(1) & d.order.eq("Chiroptera"),
            lambda d: d.order.ne("Chiroptera") | d.bat_hibernation.isna())
        print("bat matching: " + ", ".join(f"{k} {v}" for k, v in m.how.value_counts().items()))
    if mammals is not None:
        mam = sp[sp["class"].eq("Mammalia")]
        c = hib.combine(mam, hib.read_mammals(mammals), hib.read_bats(bats) if bats else None, names).set_index("slug")
        sp["mammal_hibernation"] = sp.slug.map(c.state)
        r["mammal hibernation: Heldstab et al. 2018, with Wilkinson & Adams 2019 for bats"] = (
            ("Mammalia",), lambda d: d.mammal_hibernation.eq(1), lambda d: d.mammal_hibernation.isna())
        print("mammal hibernation: " + ", ".join(f"{k} {v}" for k, v in c.how.value_counts().items()))
    tree = json.load(open(site / "src/data/generated/phylogeny.json"))["root"]
    partners = set(pd.read_parquet(site / "public/data/mitoport/partners_species.parquet").slug)
    rows = []
    for name, (classes, yes, unknown) in r.items():
        d = sp[sp["class"].isin(classes)]
        y = yes(d).fillna(False).astype(bool)
        u = unknown(d).fillna(False).astype(bool) if unknown else pd.Series(False, index=d.index)
        state = {int(i): (None if uu else int(yy)) for i, yy, uu in zip(d.id, y, u)}
        known_yes = d[y & ~u]
        rows.append({"candidate": name, "classes": ",".join(classes), "species": len(d), "yes": len(known_yes),
                     "no": int((~y & ~u).sum()), "unknown": int(u.sum()), "families_yes": known_yes.family.nunique(),
                     "trait_clades": trait_clades(tree, state), "min_changes": min_changes(tree, state),
                     "yes_with_partners": int(known_yes.slug.isin(partners).sum()),
                     "scope_with_partners": int(d.slug.isin(partners).sum())})
    t = pd.DataFrame(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(out, sep="\t", index=False)
    print(f"fishes: {fish.sum()}; FishBase match: " + ", ".join(f"{k} {v}" for k, v in match.how.value_counts().items()))
    return t


if __name__ == "__main__":
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    extra = [Path(a) for a in sys.argv[4:7]]
    print(main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), *extra).to_string(index=False))
