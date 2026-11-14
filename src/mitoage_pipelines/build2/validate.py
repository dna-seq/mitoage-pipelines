"""Checks that the species' record_status and the changelog agree (site repository: pipeline/validate.ts):
an `original` species has no changelog row, an `updated` one has at least one, a `renamed` one has a
scientific_name rename, a `retired` one a record_status row, a `new` one an added row; and every child row
of a new species is `new`."""


def problems(species: list[dict], changelog: list[dict]) -> list[str]:
    rows: dict[int, list[dict]] = {}
    for c in changelog:
        head = c["row_key"].split(":")[0]
        if head.isdigit(): rows.setdefault(int(head), []).append(c)
    out = []
    for s in species:
        r, st, sid = rows.get(s["id"], []), s["record_status"], s["id"]
        fields = {(c["table"], c["field"], c["change_type"]) for c in r}
        if st == "original" and r: out.append(f"{sid}: original but {len(r)} changelog rows")
        if st == "updated" and not r: out.append(f"{sid}: updated without a changelog row")
        if st == "renamed" and ("species", "scientific_name", "renamed") not in fields: out.append(f"{sid}: renamed without a rename row")
        if st == "retired" and ("species", "record_status", "retired") not in fields: out.append(f"{sid}: retired without a record_status row")
        if st == "new" and ("species", "(row)", "added") not in fields: out.append(f"{sid}: new without an added row")
        if st != "new" and sid >= 10_000_000: out.append(f"{sid}: a minted id on a build 1 species")
    return out


def check_children(status: dict[int, str], rows: list[dict], table: str) -> list[str]:
    return [f"{table} {r['species_id']}: row of a new species is {r['record_status']}" for r in rows
            if status.get(r["species_id"]) == "new" and r["record_status"] != "new"]
