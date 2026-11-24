"""Break each guard of `mitoport swap` in turn and check that a test of tests/test_cybrid_swap.py fails.
Prints one line per break, "caught" or "SURVIVED", and exits non-zero if any survived or was not applied.
Needs the stab environment (`uv sync --group stab --group dev`), and the inputs of the tests (structures, records, build,
vertebrate alignments), so that the tests that read them do not skip."""
import subprocess
import sys
from pathlib import Path

S = "src/mitoage_pipelines/mitoport/cybrid_swap.py"
C = "src/mitoage_pipelines/mitoport/cybrid_stab.py"
breaks = [
 (S, "the residue placed against another residue kept", "            if chain_seq[s0 + k] == human_row[cols[h0 + k]]:\n                out[s0 + k] = cols[h0 + k]", "            if True:\n                out[s0 + k] = cols[h0 + k]"),
 (S, "end gaps that cost", "extend_gap_score=-0.5, end_gap_score=0)       # a chain may be a part of the protein, at either end", "extend_gap_score=-0.5)       # a chain may be a part of the protein, at either end"),
 (S, "alignment columns taken as positions of the gapped row", "    cols = [c for c, ch in enumerate(human_row) if ch != \"-\"]", "    cols = list(range(len(human_row)))"),
 (S, "an unknown residue of the species placed", "        if col is None or aa not in THREE or row[col] not in THREE:", "        if col is None or aa not in THREE:"),
 (S, "an unknown residue of the structure placed", "        if col is None or aa not in THREE or row[col] not in THREE:", "        if col is None or row[col] not in THREE:"),
 (S, "a residue without a column read", "        if col is None or aa not in THREE or row[col] not in THREE:", "        if aa not in THREE or row[col] not in THREE:"),
 (S, "an identical residue counted as a substitution", "        if row[col] != aa:\n            subs[i] = row[col]", "        if True:\n            subs[i] = row[col]"),
 (S, "residues not placed not counted", "            lost += 1\n            continue\n        placed += 1", "            continue\n        placed += 1"),
 (S, "a residue at exactly 4.5 A not a contact", "(d[a:b] <= cs.CONTACT_A).any()", "(d[a:b] < cs.CONTACT_A).any()"),
 (S, "the query bound too short", "distance_upper_bound=cs.CONTACT_A + 1e-9", "distance_upper_bound=cs.CONTACT_A - 1e-9"),
 (S, "a contact needs every atom near", "(d[a:b] <= cs.CONTACT_A).any()", "(d[a:b] <= cs.CONTACT_A).all()"),
 (S, "only the first mitochondrial chain's atoms", "np.concatenate([chains[i][\"heavy\"] for i in mt_idx])", "chains[mt_idx[0]][\"heavy\"]"),
 (S, "a column without a residue kept", "if hit and c is not None}", "if hit}"),
 (S, "only the first copy's contact columns", "for k, f in flags.items() for c, hit", "for k, f in list(flags.items())[:1] for c, hit"),
 (S, "an unknown of the species a difference", "sum(1 for c in columns if row[c] in THREE and human_row[c] in THREE and row[c] != human_row[c])", "sum(1 for c in columns if human_row[c] in THREE and row[c] != human_row[c])"),
 (S, "a gap of the human row a difference", "sum(1 for c in columns if row[c] in THREE and human_row[c] in THREE and row[c] != human_row[c])", "sum(1 for c in columns if row[c] in THREE and row[c] != human_row[c])"),
 (S, "the host a candidate", "s not in (donor, HOST))", "s not in (donor,))"),
 (S, "the donor a candidate", "s not in (donor, HOST))", "s not in (HOST,))"),
 (S, "a non-mammal a candidate", "if s in mammals and s not in", "if s not in"),
 (S, "the window's lower bound open", "lo * n <= 4 * contact_differences", "lo * n < 4 * contact_differences"),
 (S, "the window's upper bound open", "<= hi * n], n", "< hi * n], n"),
 (S, "the candidates in the alignment's order", "    pool = sorted(s for s in rows if s in mammals", "    pool = list(s for s in rows if s in mammals"),
 (S, "the draw seeded by its number alone", "np.random.default_rng((seed, draw))", "np.random.default_rng(draw)"),
 (S, "the last candidate never drawn", ".integers(len(candidates))", ".integers(len(candidates) - 1)"),
 (S, "the sign of a term turned the wrong way", "    out[todo] = -acc / samples", "    out[todo] = acc / samples"),
 (S, "mutant minus wild type exchanged", "acc[b:b + batch] += (mut - wt_dG)", "acc[b:b + batch] += (wt_dG - mut)"),
 (S, "the generator not seeded", "    torch.manual_seed(seed)\n    acc", "    acc"),
 (S, "the wild-type pass without the shared noise", "wt_dG = model.folding_dG(domain, wt_seq, decoding_order=order, backbone_noise=noise)", "wt_dG = model.folding_dG(domain, wt_seq, decoding_order=order, backbone_noise=None)"),
 (S, "the samples not accumulated", "acc[b:b + batch] += (mut", "acc[b:b + batch] = (mut"),
 (S, "the mean over the wrong number of samples", "    out[todo] = -acc / samples", "    out[todo] = -acc / (samples + 1)"),
 (S, "binding taken as whole - one + two", "    return whole - (one + two)", "    return whole - one + two"),
 (S, "the other hosts' rows kept", "m = m[(m.host == HOST) & m[\"complex\"].isin(list(COMPLEXES)) & m.donor.isin(DONORS)]", "m = m[m[\"complex\"].isin(list(COMPLEXES)) & m.donor.isin(DONORS)]"),
 (S, "whole-cybrid rows kept", "m = m[(m.host == HOST) & m[\"complex\"].isin(list(COMPLEXES)) & m.donor.isin(DONORS)]", "m = m[(m.host == HOST) & m.donor.isin(DONORS)]"),
 (S, "donors without orthologs kept", "m = m[(m.host == HOST) & m[\"complex\"].isin(list(COMPLEXES)) & m.donor.isin(DONORS)]", "m = m[(m.host == HOST) & m[\"complex\"].isin(list(COMPLEXES))]"),
 (S, "the mitochondrial chains counted as a gene's", "        if i in mt_chain:\n            continue\n        g = chain_gene", "        g = chain_gene"),
 (S, "a chain without a gene accepted", "        if g is None:\n            raise SystemExit(f\"chain {c['id']} has no gene in chains.tsv\")\n", "        if g is None:\n            continue\n"),
 (S, "a donor without a row read", "    if donor not in rows:\n", "    if False:\n"),
 (S, "a gene without an alignment read", "        if rows is None:\n            return done(case_row(donor, complex_, gene, \"not assessed: no alignment of the gene\"))", "        if False:\n            return done(case_row(donor, complex_, gene, \"not assessed: no alignment of the gene\"))"),
 (S, "a gene without a chain taken as assessed", "    if gene not in by_gene:\n", "    if False:\n"),
 (S, "only the first copy carries the ortholog", "    u_idx, one_idx = by_gene[gene], sorted(mt_chain)", "    u_idx, one_idx = by_gene[gene][:1], sorted(mt_chain)"),
 (S, "the null run for an ortholog identical to human", "if null == \"assessed\" else []", "if len(candidates) >= MIN_CANDIDATES else []"),
 (S, "three candidates not enough", "    elif len(candidates) < MIN_CANDIDATES:", "    elif len(candidates) <= MIN_CANDIDATES:"),
 (S, "the donor's substitutions of the mitochondrial chains left out", "{**mt_sub, **extra}", "{**extra}"),
 (S, "B1 without the ortholog", "for sub in [{}, donor_sub, *species_sub.values()]:", "for sub in [{}, {}, *species_sub.values()]:"),
 (S, "restoration with the sign turned", "        rest = float(b0.mean() - b1.mean())", "        rest = float(b1.mean() - b0.mean())"),
 (S, "the null's restoration with the sign turned", "null = (b0.mean() - nul.groupby(\"draw\").ddG.mean()).to_numpy()", "null = (nul.groupby(\"draw\").ddG.mean() - b0.mean()).to_numpy()"),
 (S, "the 95th percentile replaced", "\"null_p95\": np.percentile(null, 95) if len(null) else np.nan", "\"null_p95\": np.percentile(null, 90) if len(null) else np.nan"),
 (S, "equal to the percentile exceeds it", "rest > np.percentile(null, 95))", "rest >= np.percentile(null, 95))"),
 (S, "the share of draws without the real one", "(1 + np.sum(null >= rest)) / (1 + len(null))", "(1 + np.sum(null > rest)) / (1 + len(null))"),
 (S, "ranked from the smallest", "rank(method=\"min\", ascending=False)", "rank(method=\"min\", ascending=True)"),
 (S, "ties ranked in order", "rank(method=\"min\", ascending=False)", "rank(method=\"first\", ascending=False)"),
 (S, "a restoration of zero among the largest", "(res.restoration > 0) & (res[\"rank\"]", "(res.restoration >= 0) & (res[\"rank\"]"),
 (S, "one fewer among the largest", "res[\"rank\"] <= res.genes_ranked.map(top_size)", "res[\"rank\"] < res.genes_ranked.map(top_size)"),
 (S, "a quarter not rounded up", "-(-n_genes // 4)", "n_genes // 4"),
 (S, "the limit of twelve moved", "3 if n_genes <= 12", "3 if n_genes <= 13"),
 (S, "exactly half a majority", "2 * both > len(scored)", "2 * both >= len(scored)"),
 (S, "any donor enough", "2 * both > len(scored)", "both > 0"),
 (S, "donors without an ortholog counted", "2 * both > len(scored)", "2 * both > len(c)"),
 (S, "either condition enough", "both = int((r.in_top & r.exceeds_null).sum())", "both = int((r.in_top | r.exceeds_null).sum())"),
 (S, "the union of genes summed", "common = set.intersection(*ok) if ok else set()", "common = set.union(*ok) if ok else set()"),
 (S, "a correlation of two cases", "if len(g) >= 3 and g.pct.nunique()", "if len(g) >= 2 and g.pct.nunique()"),
 (S, "a missing task not said", "    if missing:\n        print(f\"MISSING:", "    if False:\n        print(f\"MISSING:"),
 (S, "the null rows read as missing values", "scores.append(pd.read_csv(f\"{f}.tsv\", sep=\"\\t\", keep_default_na=False, na_values=[\"\"]))", "scores.append(pd.read_csv(f\"{f}.tsv\", sep=\"\\t\"))"),
 (C, "the residues' atom counts wrong", "hv += atoms; nh.append(len(atoms))", "hv += atoms; nh.append(1)"),
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
        r = subprocess.run(["uv", "run", "--frozen", "--no-sync", "pytest", "-q", "-x", "-p", "no:cacheprovider", "tests/test_cybrid_swap.py"], capture_output=True, text=True)
    finally:
        T.write_text(orig)
    ok = r.returncode != 0
    caught += ok
    print("caught " if ok else "SURVIVED", name, flush=True)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) else 1)
