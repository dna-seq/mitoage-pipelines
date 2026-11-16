"""Break each guard of H26's trait (src/mitoage_pipelines/adaptations/longevity.py) in turn and check that a test fails.

The loop is slurm/adapt-mutations.py's, with this file's mutations and tests: one line of the source replaced,
tests/test_adapt_longevity.py run, the file put back byte for byte. Exits non-zero if a mutation survived.
"""
import importlib.util
import sys
from pathlib import Path

LO = Path("src/mitoage_pipelines/adaptations/longevity.py")
spec = importlib.util.spec_from_file_location("adapt_mutations", Path(__file__).with_name("adapt-mutations.py"))
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
base.TESTS = ["tests/test_adapt_longevity.py"]
base.MUTATIONS = [
    (LO, "the middle of the residuals counts as not long-lived", "    state[has & resid.lt(mid)] = 0", "    state[has & resid.lt(hi)] = 0"),
    (LO, "a species without a residual counts as not long-lived", "    state = pd.Series(pd.NA, index=resid.index, dtype=\"Int64\")", "    state = pd.Series(0, index=resid.index, dtype=\"Int64\")"),
    (LO, "the upper cut is the median", "UPPER_Q, LOWER_Q = 0.75, 0.5", "UPPER_Q, LOWER_Q = 0.5, 0.5"),
    (LO, "a value on the upper cut is left out", "    state[has & resid.ge(hi)] = 1", "    state[has & resid.gt(hi)] = 1"),
    (LO, "species without an adult mass get a residual", "default_dataset(species) & species.adult_weight_g.gt(0) & species.slug.isin(usable)", "default_dataset(species) & species.slug.isin(usable)"),
    (LO, "lifespan records of any quality get a residual", "    m = species[\"class\"].eq(cls) & default_dataset(species) & ", "    m = species[\"class\"].eq(cls) & species.record_status.ne(\"retired\") & "),
    (LO, "species without a tip or alignments get a residual", " & species.adult_weight_g.gt(0) & species.slug.isin(usable)", " & species.adult_weight_g.gt(0)"),
    (LO, "the residual is lifespan itself", "        d[\"residual\"] = d.log10_lifespan - (intercept + slope * d.log10_mass)", "        d[\"residual\"] = d.log10_lifespan"),
    (LO, "small orders get a residual within the order", "    big = d.order.map(d.order.value_counts()).ge(MIN_ORDER) & ", "    big = d.order.map(d.order.value_counts()).ge(1) & "),
    (LO, "the residual within the order ignores mass", "        d.loc[big, \"residual_in_order\"] = cl - slope * cm", "        d.loc[big, \"residual_in_order\"] = cl"),
    (LO, "the residual within the order is the class-wide one", "        for trait, col in ((TRAIT, \"residual\"), (TRAIT_IN_ORDER, \"residual_in_order\")):", "        for trait, col in ((TRAIT, \"residual\"), (TRAIT_IN_ORDER, \"residual\")):"),
    (LO, "the other class's species enter the fit", "    m = species[\"class\"].eq(cls) & default_dataset(species)", "    m = species[\"class\"].isin(CLASSES) & default_dataset(species)"),
    (LO, "retired records get a row", "    sp = live(species)\n    parts, res = [], []", "    sp = species\n    parts, res = [], []"),
    (LO, "the selection inputs take every class", "    scope = set(species.slug[species[\"class\"].eq(cls)])\n    return selection.prepare(", "    scope = None\n    return selection.prepare("),
    (LO, "a class outside H26 is prepared", "    if cls not in CLASSES:\n        raise SystemExit(f\"--class is one of", "    if False:\n        raise SystemExit(f\"--class is one of"),
]

if __name__ == "__main__":
    sys.exit(base.main())
