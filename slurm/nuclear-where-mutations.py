"""Break each guard of `mitoport nuclear-where` in turn and check that a test of tests/test_nuclear_where.py fails.
Prints one line per break, "caught" or "SURVIVED", and exits non-zero if any survived or was not applied."""
import subprocess
import sys
from pathlib import Path

T = Path("src/mitoage_pipelines/mitoport/nuclear_where.py")
breaks = [
 ("a mitochondrial partner counted as nuclear (the residues)", "    c = contacts[contacts.partner_is_mt == 0].copy()\n    c[\"nuclear_gene\"]", "    c = contacts.copy()\n    c[\"nuclear_gene\"]"),
 ("a mitochondrial partner counted as nuclear (the masks)", "    c = contacts[contacts.partner_is_mt == 0].copy()\n    c[\"nuc\"]", "    c = contacts.copy()\n    c[\"nuc\"]"),
 ("a residue counted once per difference", "    return out.drop_duplicates([\"host\", \"donor\", \"complex\", \"nuclear_gene\", \"partner_resnum\"])", "    return out"),
 ("disagreeing papers always left out", "            if disagree == \"leave_out\":\n                continue\n", "            continue\n"),
 ("disagreeing papers kept as they are", "        if s == \"disagree\":\n            if disagree", "        if False:\n            if disagree"),
 ("a residue touching two sites counted twice", "        d[r.nuc] = d.get(r.nuc, 0) | (1 << bit[r.nuc][r.partner_resnum])", "        d[r.nuc] = d.get(r.nuc, 0) + (1 << bit[r.nuc][r.partner_resnum])"),
 ("the union of residues over sites taken as the last site's", "            acc[u] = acc.get(u, 0) | m", "            acc[u] = m"),
 ("the area with ties counted zero", "((a > b) + 0.5 * (a == b)).mean(axis=(0, 1))", "(a > b).mean(axis=(0, 1))"),
 ("reduced and not reduced exchanged in the area", "    a, b = counts[red][:, None, :], counts[~red][None, :, :]", "    a, b = counts[~red][:, None, :], counts[red][None, :, :]"),
 ("the drawn sites not kept per protein", "        out[k] = [(g, int(p)) for g, m in n.items() for p in rng.choice(universe[g], size=m, replace=False)]", "        out[k] = [(g, int(p)) for g, m in n.items() for p in rng.choice(sorted({q for u in universe.values() for q in u}), size=m, replace=False)]"),
 ("sites drawn with replacement", "size=m, replace=False)]", "size=m, replace=True)]"),
 ("a different seed for each run", "        rng = np.random.default_rng(seed)\n        null", "        rng = np.random.default_rng()\n        null"),
 ("the subunit p-value without the observed draw", "(1 + np.sum(null[:, k] >= auc[k])) / (1 + draws) if not np.isnan(auc[k]) else np.nan", "np.sum(null[:, k] >= auc[k]) / draws if not np.isnan(auc[k]) else np.nan"),
 ("the family-wise p over the subunit's own area", "np.nanmax(null, axis=1) >= best", "np.nanmax(null, axis=1) >= np.nanmin(auc)"),
 ("a subunit named without the family-wise rule", "\"named\": bool(fw < 0.05 and auc[k] == best)", "\"named\": bool(auc[k] == best)"),
 ("donors not left out", "        lodo = np.array([auc_columns(obs[donors != d], red[donors != d]) for d in sorted(set(donors))])", "        lodo = np.array([auc_columns(obs, red) for d in sorted(set(donors))])"),
 ("every subunit in one profile group", "            same.setdefault(tuple(obs[:, k]), []).append(u)", "            same.setdefault(0, []).append(u)"),
 ("a pair without a difference dropped from the table", "[[pair_counts(ss[k], masks).get(u, 0) for u in subs] for k in pairs]", "[[pair_counts(ss[k], masks).get(u, 0) for u in subs] for k in pairs if ss[k]]"),
 ("the share of reduced over the wrong total", "g[\"share_reduced\"] = [r.reduced / n_red.get(r.complex, np.nan) for r in g.itertuples()]", "g[\"share_reduced\"] = [r.reduced / n_not.get(r.complex, np.nan) for r in g.itertuples()]"),
 ("the ranking not by the difference", "    return g.sort_values([\"difference\", \"reduced\"], ascending=False).reset_index(drop=True)", "    return g.sort_values([\"reduced\"], ascending=False).reset_index(drop=True)"),
]
caught = 0
for name, a, b in breaks:
    orig = T.read_text()
    if orig.count(a) != 1:
        print("NOT APPLIED", name); continue
    try:
        T.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "pytest", "-q", "-x", "-p", "no:cacheprovider", "tests/test_nuclear_where.py"], capture_output=True, text=True)
    finally:
        T.write_text(orig)
    ok = r.returncode != 0
    caught += ok
    print("caught " if ok else "SURVIVED", name)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) else 1)
