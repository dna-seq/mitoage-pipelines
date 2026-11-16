"""Hibernation sources and their matching to the build (adaptations/hibernation.py), on the files of
`download hibernation`, the downloaded build and NCBI's names.dmp (MITOAGE_TAXDUMP: the folder holding it)."""
import os
import shutil
from pathlib import Path

import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.adaptations import hibernation as hib
from mitoage_pipelines.build2 import pins

BATS = paths.HIBERNATION / hib.BATS_FILE
NAMES = Path(os.environ.get("MITOAGE_TAXDUMP", paths.B2_TAXDUMP / pins.TAXDUMP_BUILD2)) / "names.dmp"
need_bats = pytest.mark.skipif(not BATS.is_file(), reason="needs `download hibernation`")
need_all = pytest.mark.skipif(not (BATS.is_file() and NAMES.is_file() and (paths.BUILD / "species.parquet").is_file()),
                              reason="needs `download hibernation`, `download build` and NCBI's names.dmp")


def lines(path: Path) -> list[str]:
    text = path.read_bytes().decode().replace("\r\n", "\n").replace("\r", "\n")
    return [x for x in text.split("\n") if x.strip()]


@pytest.fixture(scope="module")
def bats():
    return hib.read_bats(BATS)


@pytest.fixture(scope="module")
def build_mammals():
    sp = pd.read_parquet(paths.BUILD / "species.parquet")
    return sp[sp.record_status.ne("retired") & sp["class"].eq("Mammalia")].reset_index(drop=True)


@need_bats
def test_bats_every_row_read(bats):
    raw = lines(BATS)
    head = raw[0].split("\t")
    col = head.index("Hibernation")
    values = [r.split("\t")[col].strip() for r in raw[1:]]
    assert len(bats) == len(raw) - 1
    assert bats.hibernation.sum() == values.count("Yes")
    assert (bats.hibernation == 0).sum() == values.count("No")
    assert bats.name.is_unique and bats.name.str.contains(" ").all()


@need_bats
@pytest.mark.parametrize("change", ["value", "duplicate", "empty"])
def test_bats_refused(tmp_path, change):
    raw = lines(BATS)
    col = raw[0].split("\t").index("Hibernation")
    first = raw[1].split("\t")
    if change == "value":
        first[col] = "Sometimes"
        raw[1] = "\t".join(first)
    elif change == "duplicate":
        raw.append(raw[1])
    else:
        first[raw[0].split("\t").index("Genus_species")] = " "
        raw[1] = "\t".join(first)
    p = tmp_path / "bats.txt"
    p.write_text("\r".join(raw))
    with pytest.raises(hib.Refused):
        hib.read_bats(p)


@need_all
def test_match_bats(bats, build_mammals):
    sp = build_mammals[build_mammals.order.eq("Chiroptera")]
    m = hib.match(sp, bats, NAMES).set_index("slug")
    table = dict(zip(bats.name, bats.hibernation))
    exact = sp[sp.scientific_name.isin(table)]
    assert set(m.index[m.how == "exact name"]) == set(exact.slug)
    for slug, sci in zip(exact.slug, exact.scientific_name):
        assert m.at[slug, "state"] == table[sci]
    # a synonym match is a table name NCBI gives to the species' own taxid, with that row's state
    tids = hib.ncbi_taxids(NAMES, set(bats.name))
    for slug, row in m[m.how == "NCBI synonym"].iterrows():
        tid = int(sp.set_index("slug").at[slug, "ncbi_taxid"])
        for name in row.table_name.split("; "):
            assert tid in tids[name] and table[name] == row.state
    # every bat of the build without such a name is unknown, not 0
    reach = {t for n, ts in tids.items() for t in ts}
    lost = sp[~sp.scientific_name.isin(table) & ~sp.ncbi_taxid.isin(reach)]
    assert set(m.index[m.state.isna()]) == set(lost.slug)


@need_all
def test_absent_is_unknown_not_no(bats, build_mammals):
    m = hib.match(build_mammals, bats, NAMES)
    assert m.state.isna().sum() == (m.how != "exact name").sum() - (m.how == "NCBI synonym").sum()
    assert m.state.notna().sum() <= len(bats)
    # a 0 comes only from a "No" in the table
    no = set(bats.name[bats.hibernation == 0])
    zero = m[m.state == 0]
    assert len(zero) and all(set(t.split("; ")) <= no or sci in no for t, sci in
                             zip(zero.table_name, build_mammals.set_index("slug").loc[zero.slug, "scientific_name"]))


@need_all
def test_ambiguous_when_rows_disagree(bats, build_mammals):
    """Two table names that NCBI gives to one build species, with opposite states, leave it unknown."""
    sp = build_mammals[build_mammals.order.eq("Chiroptera")]
    m = hib.match(sp, bats, NAMES)
    syn = m[m.how == "NCBI synonym"].iloc[0]
    tid = int(sp.set_index("slug").at[syn.slug, "ncbi_taxid"])
    other = None
    with open(NAMES, encoding="utf-8") as f:
        for line in f:
            p = line.split("\t|\t")
            if int(p[0]) == tid and p[1] not in set(bats.name) and p[1] not in set(sp.scientific_name):
                other = p[1]
                break
    if other is None:
        pytest.skip("NCBI gives this taxon no third name")
    extra = pd.DataFrame([{"name": other, "family": "", "hibernation": 1 - int(syn.state)}])
    m2 = hib.match(sp, pd.concat([bats, extra], ignore_index=True), NAMES).set_index("slug")
    assert pd.isna(m2.at[syn.slug, "state"]) and m2.at[syn.slug, "how"].startswith("ambiguous")


@need_bats
def test_fetch_checks_both_files(tmp_path, monkeypatch):
    from mitoage_pipelines.build2 import download as dl
    monkeypatch.setattr(dl, "fetch_to", lambda url, dest: shutil.copy(BATS, dest))
    manifest = tmp_path / "MANIFEST.json"
    got = hib.fetch(tmp_path / "h", manifest)
    assert got["mammals"] is None and got["bats"].read_bytes() == BATS.read_bytes()
    (tmp_path / "h" / hib.MAMMALS_FILE).write_bytes(b"not the Dryad file")
    with pytest.raises(hib.Refused):
        hib.fetch(tmp_path / "h", manifest)
    monkeypatch.setattr(dl, "fetch_to", lambda url, dest: Path(dest).write_bytes(BATS.read_bytes() + b"\r"))
    with pytest.raises(hib.Refused):
        hib.fetch(tmp_path / "h2", manifest)


@need_all
def test_ambiguous_when_one_name_reaches_two_species(bats, build_mammals):
    """A table name reached through one taxid by two build rows (here a copy of a bat under another name) leaves
    both unknown."""
    sp = build_mammals[build_mammals.order.eq("Chiroptera")]
    syn = hib.match(sp, bats, NAMES).query("how == 'NCBI synonym'").iloc[0]
    row = sp[sp.slug == syn.slug].copy()
    row["slug"], row["scientific_name"] = row.slug + "-copy", "Copyus duplicatus"
    m = hib.match(pd.concat([sp, row], ignore_index=True), bats, NAMES).set_index("slug")
    for s in (syn.slug, syn.slug + "-copy"):
        assert pd.isna(m.at[s, "state"]) and m.at[s, "how"].startswith("ambiguous")
