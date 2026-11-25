"""Break each guard of `mitoport saprot` in turn and check that a test of tests/test_cybrid_saprot.py (and, for the swap's code it reads,
tests/test_cybrid_swap.py) fails. Prints one line per break, "caught" or "SURVIVED", and exits non-zero if any survived or was not applied.
Needs the saprot and stab environments (`uv sync --group saprot --group stab --group dev`) and the inputs of the tests (structures, records,
build, vertebrate alignments, the structure tokens), so that the tests that read them do not skip."""
import subprocess
import sys
from pathlib import Path

S = "src/mitoage_pipelines/mitoport/cybrid_saprot.py"
breaks = [
 # the structure tokens
 (S, "the states kept in upper case", "out[chain] = (f[1], f[2].lower())", "out[chain] = (f[1], f[2])"),
 (S, "a line of two fields accepted", "        if len(f) < 3:\n            raise SystemExit(f\"a line of Foldseek's output has", "        if False:\n            raise SystemExit(f\"a line of Foldseek's output has"),
 (S, "a chain seen twice accepted", "        if chain in out:\n            raise SystemExit", "        if False:\n            raise SystemExit"),
 (S, "residues and states of different lengths accepted", "        if len(f[1]) != len(f[2]):\n            raise SystemExit", "        if False:\n            raise SystemExit"),
 (S, "the chain named by the text before the underscore", "chain = f[0].split(\" \")[0].rsplit(\"_\", 1)[1]", "chain = f[0].split(\" \")[0].rsplit(\"_\", 1)[0]"),
 (S, "the states' length not checked", "if len(tok_aa) != len(chain_seq) or len(tok_di) != len(chain_seq):", "if len(tok_aa) != len(chain_seq):"),
 (S, "the residues' length not checked", "if len(tok_aa) != len(chain_seq) or len(tok_di) != len(chain_seq):", "if len(tok_di) != len(chain_seq):"),
 (S, "a different residue accepted", "if any(a != b and b != \"X\" for a, b in zip(chain_seq, tok_aa)):", "if False:"),
 (S, "Foldseek's X not accepted", "if any(a != b and b != \"X\" for a, b in zip(chain_seq, tok_aa)):", "if any(a != b for a, b in zip(chain_seq, tok_aa)):"),
 (S, "the chain's own X accepted against a residue", "if any(a != b and b != \"X\" for a, b in zip(chain_seq, tok_aa)):", "if any(a != b and b != \"X\" and a != \"X\" for a, b in zip(chain_seq, tok_aa)):"),
 (S, "no binary and no cache accepted", "    if binary is None or not binary.is_file():\n        raise SystemExit", "    if False:\n        raise SystemExit"),
 (S, "the chain table not matched by id", "return {i: chain_states(c[\"seq\"], *tokens[c[\"id\"]]) if c[\"id\"] in tokens else None for i, c in enumerate(chains)}", "return {i: chain_states(c[\"seq\"], *list(tokens.values())[0]) for i, c in enumerate(chains)}"),
 # the input
 (S, "a foreign vocabulary accepted", "if len(set(v)) != len(v) or v[:5] != [\"<cls>\", \"<pad>\", \"<eos>\", \"<unk>\", \"<mask>\"]:", "if False:"),
 (S, "a residue that is no amino acid read as itself", "return [(a if a in AA else \"#\") + s for a, s in zip(seq, states)]", "return [a + s for a, s in zip(seq, states)]"),
 (S, "residues and states of different lengths read", "    if len(seq) != len(states):\n        raise ValueError", "    if False:\n        raise ValueError"),
 (S, "a place without an amino acid substituted", "if s[i] not in AA or aa not in AA:", "if aa not in AA:"),
 (S, "a replacement that is no amino acid accepted", "if s[i] not in AA or aa not in AA:", "if s[i] not in AA:"),
 (S, "the linker one residue short", "[\"G#\"] * LINKER + u, len(m) + LINKER", "[\"G#\"] * (LINKER - 1) + u, len(m) + LINKER - 1"),
 (S, "the linker given a structure", "[\"G#\"] * LINKER + u, len(m) + LINKER", "[\"Gp\"] * LINKER + u, len(m) + LINKER"),
 (S, "the offset without the linker", "[\"G#\"] * LINKER + u, len(m) + LINKER", "[\"G#\"] * LINKER + u, len(m)"),
 (S, "the chains in the other order", "return m + [\"G#\"] * LINKER + u, len(m) + LINKER", "return u + [\"G#\"] * LINKER + m, len(m) + LINKER"),
 (S, "the mitochondrial substitutions put on the nuclear chain", "m = sa_tokens(apply(m_seq, m_subs or {}), m_states)\n    u = sa_tokens(apply(u_seq, u_subs or {}), u_states)", "m = sa_tokens(apply(m_seq, u_subs or {}), m_states)\n    u = sa_tokens(apply(u_seq, m_subs or {}), u_states)"),
 (S, "the shuffle seeded by the draw alone", "np.random.default_rng((1, draw)).permutation", "np.random.default_rng((0, draw)).permutation"),
 (S, "the shuffle over every place", "idx = [i for i, a in enumerate(u_seq) if a in AA]", "idx = list(range(len(u_seq)))"),
 (S, "the shuffle listing unchanged places", "for k, p in enumerate(perm) if u_seq[idx[p]] != u_seq[idx[k]]}", "for k, p in enumerate(perm)}"),
 # the pseudo-log-likelihood
 (S, "the masked place left as it was", "x[np.arange(len(part)), cols] = [vocab[\"#\" + s] for s in state]", "x[np.arange(len(part)), cols] = x[np.arange(len(part)), cols]"),
 (S, "the masked place without its structure state", "x[np.arange(len(part)), cols] = [vocab[\"#\" + s] for s in state]", "x[np.arange(len(part)), cols] = [vocab[\"##\"] for s in state]"),
 (S, "the model's own mask token", "x[np.arange(len(part)), cols] = [vocab[\"#\" + s] for s in state]", "x[np.arange(len(part)), cols] = [vocab[\"<mask>\"] for s in state]"),
 (S, "the column without the start token", "cols = np.array([positions[k] + 1 for k in part])", "cols = np.array([positions[k] for k in part])"),
 (S, "no end token", "ids = np.array([vocab[\"<cls>\"], *(vocab[t] for t in tokens), vocab[\"<eos>\"]], dtype=np.int64)", "ids = np.array([vocab[\"<cls>\"], *(vocab[t] for t in tokens)], dtype=np.int64)"),
 (S, "not normalised over the twenty", "out[k] = (lp - logsumexp(lp))[AA.index(tokens[positions[k]][0])]", "out[k] = lp[AA.index(tokens[positions[k]][0])]"),
 (S, "the wrong residue's probability", "out[k] = (lp - logsumexp(lp))[AA.index(tokens[positions[k]][0])]", "out[k] = (lp - logsumexp(lp))[(AA.index(tokens[positions[k]][0]) + 1) % 20]"),
 (S, "the twenty at another state", "lp = logits[j, [vocab[a + state[j]] for a in AA]]", "lp = logits[j, [vocab[a + state[(j + 1) % len(state)]] for a in AA]]"),
 (S, "a place without an amino acid read", "todo = [k for k, p in enumerate(positions) if tokens[p][0] in AA]", "todo = list(range(len(positions)))"),
 (S, "only the first batch read", "    for b in range(0, len(todo), batch):", "    for b in range(0, min(len(todo), batch), batch):"),
 (S, "the batch's rows mixed up", "        for j, k in enumerate(part):\n            lp", "        for j, k in enumerate(reversed(part)):\n            lp"),
 # the pairs
 (S, "a contact needs the first residue only", "            out[k] = [i for i, f in enumerate(flags[k]) if f]", "            out[k] = [i for i, f in enumerate(flags[k]) if f][:1]"),
 (S, "the pairs of chains that touch none kept", "        if flags.any():\n            out[k]", "        if True:\n            out[k]"),
 (S, "the distance to the first atom only", "d[a:b].min() for a, b in zip(edges[:-1], edges[1:])", "d[a] for a, b in zip(edges[:-1], edges[1:])"),
 (S, "the distance to the first mitochondrial chain only", "cKDTree(np.concatenate([chains[i][\"heavy\"] for i in mt_idx]))", "cKDTree(chains[mt_idx[0]][\"heavy\"])"),
 # one task
 (S, "a task with a table read again", "    if Path(f\"{stem}.cases.tsv\").is_file():\n        print(\"skipped", "    if False:\n        print(\"skipped"),
 (S, "the marker written first", "        pd.DataFrame([case]).to_csv(f\"{stem}.cases.tmp\"", "        raise SystemExit('x')\n        pd.DataFrame([case]).to_csv(f\"{stem}.cases.tmp\""),
 (S, "a table written in place", "            Path(f\"{stem}{name}.tmp\").replace(f\"{stem}{name}.tsv\")", "            Path(f\"{stem}{name}.tmp\").rename(f\"{stem}{name}.tsv\"); Path(f\"{stem}{name}.tmp\").write_text('')"),
 (S, "the pair longer than the model read", "if len(chains[k][\"seq\"]) + LINKER + len(chains[c][\"seq\"]) > MAX_RESIDUES:", "if False:"),
 (S, "a pair of exactly the length refused", "if len(chains[k][\"seq\"]) + LINKER + len(chains[c][\"seq\"]) > MAX_RESIDUES:", "if len(chains[k][\"seq\"]) + LINKER + len(chains[c][\"seq\"]) >= MAX_RESIDUES:"),
 (S, "the linker not counted in the length", "if len(chains[k][\"seq\"]) + LINKER + len(chains[c][\"seq\"]) > MAX_RESIDUES:", "if len(chains[k][\"seq\"]) + len(chains[c][\"seq\"]) > MAX_RESIDUES:"),
 (S, "a chain without states read", "            elif states[k] is None or states[c] is None:", "            elif False:"),
 (S, "only the nuclear chain's states checked", "            elif states[k] is None or states[c] is None:", "            elif states[c] is None:"),
 (S, "no pair said as not assessed", "    if not pairs:\n        why", "    if False:\n        why"),
 (S, "the null drawn for a donor that changes nothing", "drawn = [sw.draw_species(candidates, k) for k in range(draws)] if null == \"assessed\" else []", "drawn = [sw.draw_species(candidates, k) for k in range(draws)] if len(candidates) >= sw.MIN_CANDIDATES else []"),
 (S, "the shuffles of a donor that changes nothing read", "for j in range(shuffles)} if donor_sub else {}", "for j in range(shuffles)}"),
 (S, "B1 read without the donor's ortholog", "for sub in [{}, donor_sub, *species_sub.values(), *shuf.values(), *single.values()]:", "for sub in [{}, {}, *species_sub.values(), *shuf.values(), *single.values()]:"),
 (S, "B1 labelled with the donor's substitutions left out", "res = [row(\"B0\", {}, -1, host), row(\"B1\", donor_sub, -1, donor),", "res = [row(\"B0\", {}, -1, host), row(\"B1\", {}, -1, donor),"),
 (S, "the donor's substitutions left off the mitochondrial chain", "mt_by_chain = {k: {i: aa for (kk, i), aa in mt_sub.items() if kk == k} for k, _, _ in pairs}", "mt_by_chain = {k: {} for k, _, _ in pairs}"),
 (S, "the nuclear chain's substitutions put on every chain", "{i: aa for (cc, i), aa in sub.items() if cc == c}) for k, c, pos in pairs}", "{i: aa for (cc, i), aa in sub.items()}) for k, c, pos in pairs}"),
 (S, "the misfit with the sign of the PLL", "\"pll\": pll, \"ddG\": -pll}", "\"pll\": pll, \"ddG\": pll}"),
 (S, "the PLL of the first pair only", "total = lambda sub: float(sum(np.nansum(v) for v in logp[key(sub)].values()))", "total = lambda sub: float(sum(np.nansum(v) for v in list(logp[key(sub)].values())[:1]))"),
 (S, "the positions counted from B1", "n_pos = int(sum(np.isfinite(v).sum() for v in logp[key({})].values()))", "n_pos = 0"),
 (S, "a residue's effect against the donor's whole chain", "delta = logp[key(single[(c, i)])][(k, c)] - logp[key({})][(k, c)]", "delta = logp[key(single[(c, i)])][(k, c)] - logp[key(donor_sub)][(k, c)]"),
 (S, "a residue's effect with the sign turned", "delta = logp[key(single[(c, i)])][(k, c)] - logp[key({})][(k, c)]", "delta = logp[key({})][(k, c)] - logp[key(single[(c, i)])][(k, c)]"),
 (S, "a residue read beside every chain of the gene", "            if cc != c:\n                continue\n            delta", "            delta"),
 (S, "the interface by 4.6 A", "\"interface\": bool(dist[c][i] <= cs.CONTACT_A)", "\"interface\": bool(dist[c][i] <= cs.CONTACT_A + 0.1)"),
 (S, "the interface by less than 4.5 A", "\"interface\": bool(dist[c][i] <= cs.CONTACT_A)", "\"interface\": bool(dist[c][i] < cs.CONTACT_A)"),
 (S, "a residue's number taken from its index", "\"resnum\": chains[c][\"num\"][i]", "\"resnum\": i + 1"),
 (S, "the host's own subunit read as another host's", "if host != sw.HOST and donor_source is None:\n        raise SystemExit", "if False:\n        raise SystemExit"),
 # the tables
 (S, "the real partner against the shuffled the wrong way", "\"real_minus_shuffled\": b0 - sh.pll.mean() if len(sh) else np.nan", "\"real_minus_shuffled\": sh.pll.mean() - b0 if len(sh) else np.nan"),
 (S, "a residue's effect not summed over its chains", "g.delta_pll.sum()], axis=1)", "g.delta_pll.first()], axis=1)"),
 (S, "the residues ranked from the smallest", "t[i].groupby([\"donor\", \"complex\", \"gene\"]).delta_pll.rank(method=\"min\", ascending=False)", "t[i].groupby([\"donor\", \"complex\", \"gene\"]).delta_pll.rank(method=\"min\", ascending=True)"),
 (S, "the residues away from the interface ranked", "t.loc[i, \"rank_in_gene\"] = t[i].groupby", "t[\"rank_in_gene\"] = t.groupby"),
 (S, "the signed change compared at the interface", "a = t.delta_pll.abs()", "a = t.delta_pll"),
 (S, "the interface compared the wrong way", "u = stats.mannwhitneyu(x, y, alternative=\"greater\")", "u = stats.mannwhitneyu(x, y, alternative=\"less\")"),
 (S, "the pairs formed across complexes", "    for _, g in df.groupby([\"host\", \"complex\"]):", "    for _, g in [(0, df)]:"),
 (S, "ties counted as wins", "wins += float(((x > nr) + 0.5 * (x == nr)).sum())", "wins += float(((x >= nr)).sum())"),
 (S, "ties counted as losses", "wins += float(((x > nr) + 0.5 * (x == nr)).sum())", "wins += float(((x > nr)).sum())"),
 (S, "the pairs counted twice", "wins += float(((x > nr) + 0.5 * (x == nr)).sum()); n += len(nr)", "wins += float(((x > nr) + 0.5 * (x == nr)).sum()); n += 2 * len(nr)"),
 (S, "the genes of one donor summed", "common = set.intersection(*(set(g[g.donor == d].gene) & set(stb[(stb.donor == d) & (stb[\"complex\"] == cx)].gene) for d in donors))", "common = set(g[g.donor == donors[0]].gene) & set(stb[(stb.donor == donors[0]) & (stb[\"complex\"] == cx)].gene)"),
 (S, "StaB-ddG's genes not restricted to the common ones", "(stb[\"complex\"] == cx) & stb.gene.isin(common)].restoration.sum()", "(stb[\"complex\"] == cx)].restoration.sum()"),
 (S, "the count over every subunit", "float(e[have].sum(axis=1).iloc[0]) if len(e) and have else np.nan", "float(e.drop(columns=[\"host\", \"donor\", \"complex\", \"status\"]).sum(axis=1).iloc[0]) if len(e) and have else np.nan"),
 (S, "the other host's step 1 rows read", "e = expo[(expo.host == host) & (expo.donor == d) & (expo[\"complex\"] == cx)]", "e = expo[(expo.donor == d) & (expo[\"complex\"] == cx)]"),
 (S, "a missing task not said", "    if missing:\n        print(f\"MISSING:", "    if False:\n        print(f\"MISSING:"),
 (S, "the null rows read as missing values", "pd.read_csv(f\"{f}.{name}.tsv\" if name else f\"{f}.tsv\", sep=\"\\t\", keep_default_na=False, na_values=[\"\"])", "pd.read_csv(f\"{f}.{name}.tsv\" if name else f\"{f}.tsv\", sep=\"\\t\")"),
 (S, "a missing task read as zero", "        if not Path(f\"{f}.cases.tsv\").is_file():\n            missing.append", "        if False:\n            missing.append"),
]
only = sys.argv[1:]                       # names (or parts of names) of the breaks to run, all when none is given
breaks = [b for b in breaks if not only or any(o in b[1] for o in only)]
caught = 0
for name_file, name, a, b in breaks:
    T = Path(name_file)
    orig = T.read_text()
    if orig.count(a) != 1:
        print("NOT APPLIED", name); continue
    try:
        T.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "--no-sync", "pytest", "-q", "-x", "-p", "no:cacheprovider", "tests/test_cybrid_saprot.py"], capture_output=True, text=True)
    finally:
        T.write_text(orig)
    ok = r.returncode != 0
    caught += ok
    print("caught " if ok else "SURVIVED", name, flush=True)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) else 1)
