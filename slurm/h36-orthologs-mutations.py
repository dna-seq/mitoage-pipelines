"""Break each guard of the donor orthologs (`mitoport/donor_orthologs.py`) and of the swap's reading of a donor's protein
(`cybrid_swap.py`: `place_protein`, `protein_subs`, `protein_contact_differences`, the `donor_source` branch of `run`) in turn and
check that a test fails. Prints one line per break, "caught" or "SURVIVED", and exits non-zero if any survived or was not applied.
The donor tests need only the dev group; the swap's need the stab environment and the inputs of the swap's tests (structures,
records, build, vertebrate alignments, `mitoage download donor-orthologs`), so that none skips:
  UV_PROJECT_ENVIRONMENT=.venv-stab uv run --frozen --no-sync python slurm/h36-orthologs-mutations.py"""
import subprocess
import sys
from pathlib import Path

D = "src/mitoage_pipelines/mitoport/donor_orthologs.py"
S = "src/mitoage_pipelines/mitoport/cybrid_swap.py"
TD = ["tests/test_donor_orthologs.py"]
TS = ["tests/test_cybrid_swap.py", "-k", "ProteinReading or ProteinTask"]
OLD_REF = 'r.get("assembly_info", {}).get("refseq_category") == "reference genome" and r["accession"].startswith("GCF_")'
CHECK = 'if len(s) != int(r.length) or sequence_digest(s) != r.sha256:'
breaks = [
 (D, TD, "a RefSeq assembly that is not the reference taken", OLD_REF, 'r["accession"].startswith("GCF_")'),
 (D, TD, "a GenBank assembly taken", OLD_REF, 'r.get("assembly_info", {}).get("refseq_category") == "reference genome"'),
 (D, TD, "two reference genomes not refused", "    if len(ref) > 1:", "    if len(ref) > 5:"),
 (D, TD, "the human gene found by a part of the symbol", 'if r["gene"]["symbol"] == symbol]', 'if symbol in r["gene"]["symbol"]]'),
 (D, TD, "the group not split by taxon", 'out.setdefault(int(r["gene"]["tax_id"]), [])', "out.setdefault(0, [])"),
 (D, TD, "the genes of a species not in gene id order", 'return {t: sorted(g, key=lambda x: int(x["gene_id"])) for t, g in out.items()}', "return {t: list(g) for t, g in out.items()}"),
 (D, TD, "the annotations of every assembly counted", 'if a.get("assembly_accession") == assembly\n', "if True\n"),
 (D, TD, "a protein on no genomic accession of the assembly taken", "        if on:\n            seen.setdefault", "        if True:\n            seen.setdefault"),
 (D, TD, "XP_ before NP_", "(0 if prefix == \"NP\" else 1, prefix, int(num), int(ver))", "(1 if prefix == \"NP\" else 0, prefix, int(num), int(ver))"),
 (D, TD, "the accession number as text", "prefix, int(num), int(ver))", "prefix, num, int(ver))"),
 (D, TD, "the version ignored", "prefix, int(num), int(ver))", "prefix, int(num), 0)"),
 (D, TD, "the shortest protein taken", "top = max(p[\"length\"] for p in proteins)", "top = min(p[\"length\"] for p in proteins)"),
 (D, TD, "a protein of equal length taken by the order listed", 'key=lambda p: accession_key(p["accession"]))', "key=lambda p: 0)"),
 (D, TD, "the proteins not taken not reported", 'return tied[0], [p["accession"] for p in tied[1:]]', "return tied[0], []"),
 (D, TD, "no protein to choose from accepted", '    if not proteins:\n        raise ValueError("no protein to choose from")', "    if False:\n        raise ValueError(\"no protein to choose from\")"),
 (D, TD, "paralogs: the first taken", "    if len(genes_of_species) > 1:", "    if len(genes_of_species) > 5:"),
 (D, TD, "a missing ortholog not said", 'return None, "no ortholog"', 'return None, "ok"'),
 (D, TD, "a gene that is not protein coding taken", '    if g.get("type") != "PROTEIN_CODING":', "    if False:"),
 (D, TD, "a gene not annotated on the assembly taken", "    if not assembly_accessions(g, assembly):", "    if False:"),
 (D, TD, "the candidates of an ambiguous gene not listed", "        if len(cand) > 1:", "        if len(cand) > 5:"),
 (D, TD, "a gene without a protein on the assembly taken", 'row["status"] = "ortholog has no RefSeq protein on the assembly"', "row[\"status\"] = OK"),
 (D, TD, "the number of proteins wrong", "proteins=len(prots)", "proteins=1"),
 (D, TD, "the length NCBI gave not checked", '        if len(s) != int(t.at[i, "length"]):', "        if False:"),
 (D, TD, "the digest of another string", 't.at[i, "sha256"] = sequence_digest(s)', 't.at[i, "sha256"] = sequence_digest(s.lower())'),
 (D, TD, "equal isoforms counted as different", "sum(seqs[a] != s for a in", "sum(seqs[a] == s for a in"),
 (D, TD, "a header twice accepted", '            if name in out:\n                raise ValueError(f"{name} twice in the FASTA")', "            if False:\n                raise ValueError(f\"{name} twice in the FASTA\")"),
 (D, TD, "only the last line of a sequence kept", "out[name] += line.strip().upper()", "out[name] = line.strip().upper()"),
 (D, TD, "the case of a sequence kept", "out[name] += line.strip().upper()", "out[name] += line.strip()"),
 (D, TD, "the length of a fetched sequence not checked", CHECK, "if sequence_digest(s) != r.sha256 and len(s) > 10 ** 6:"),
 (D, TD, "the digest of a fetched sequence not checked", CHECK, "if len(s) != int(r.length):"),
 (D, TD, "the gaps fetched", "seqs = fetch(list(t.protein))", "seqs = fetch(list(read_table(table_path).protein))"),
 (D, TD, "the files' digests left out of the manifest", "sha256=digests, table_sha256=", "sha256={}, table_sha256="),
 (D, TD, "the table's digest left out of the manifest", 'table_sha256=hashlib.sha256(table_path.read_bytes()).hexdigest())', 'table_sha256="")'),
 (D, TD, "a changed file read back", '            if sequence_digest(seqs[r.protein]) != r.sha256:\n                raise SystemExit(f"{f}: {r.protein} is not the pinned sequence")', "            if False:\n                raise SystemExit(f\"{f}: {r.protein} is not the pinned sequence\")"),
 (D, TD, "a record the table does not hold accepted", "if set(seqs) != set(g.protein):", "if not set(g.protein) <= set(seqs):"),
 (D, TD, "a missing file not refused", '        if not f.is_file():\n            raise SystemExit(f"{f} is missing', '        if False:\n            raise SystemExit(f"{f} is missing'),
 (D, TD, "an uncovered species read from the table", "        if species not in self.species:\n            return None, None", "        if False:\n            return None, None"),
 (D, TD, "a species without a source read from the table", "        if not self.species[species]:", "        if False:"),
 (D, TD, "the reason of a gap lost", 'return None, f"not assessed: no ortholog of the donor in RefSeq ({st})"', 'return None, "not assessed: no ortholog of the donor in RefSeq ()"'),
 (D, TD, "a gene missing from the table taken as a gap", "        if st is None:\n            raise SystemExit(", "        if False:\n            raise SystemExit("),
 (D, TD, "the table not in species and gene order", '.sort_values(["_s", "_g"])', '.sort_values(["_g", "_s"])'),
 # the swap's reading of a donor's protein
 (S, TS, "end gaps that cost, in the placement", "extend_gap_score=-0.5, end_gap_score=0)\n    a = al.align(chain_seq, protein.upper())[0]", "extend_gap_score=-0.5)\n    a = al.align(chain_seq, protein.upper())[0]"),
 (S, TS, "a residue set against a gap read as a residue", 'out = ["-"] * len(chain_seq)\n    for (s0, s1), (p0, p1)', 'out = ["A"] * len(chain_seq)\n    for (s0, s1), (p0, p1)'),
 (S, TS, "the donor's letters read as written", "out[s0 + k] = protein[p0 + k].upper()", "out[s0 + k] = protein[p0 + k]"),
 (S, TS, "the placement shifted by one", "out[s0 + k] = protein[p0 + k].upper()", "out[s0 + k] = protein[p0 + k - 1].upper()"),
 (S, TS, "the chain's positions not the row's", "return (*ortholog_subs(chain_seq, list(range(len(chain_seq))), row), row)", "return (*ortholog_subs(chain_seq, [None] * len(chain_seq), row), row)"),
 (S, TS, "a gap in the donor a difference", "if hit and ids[k][i] is not None and rows[k][i] in THREE and seqs[k][i] in THREE and rows[k][i] != seqs[k][i]:", "if hit and ids[k][i] is not None and seqs[k][i] in THREE and rows[k][i] != seqs[k][i]:"),
 (S, TS, "an unknown of the chain a difference", "if hit and ids[k][i] is not None and rows[k][i] in THREE and seqs[k][i] in THREE and rows[k][i] != seqs[k][i]:", "if hit and ids[k][i] is not None and rows[k][i] in THREE and rows[k][i] != seqs[k][i]:"),
 (S, TS, "a position without an identity counted", "if hit and ids[k][i] is not None and rows[k][i] in THREE", "if hit and rows[k][i] in THREE"),
 (S, TS, "a difference away from the contacts counted", "if hit and ids[k][i] is not None and rows[k][i] in THREE", "if ids[k][i] is not None and rows[k][i] in THREE"),
 (S, TS, "a position counted once per copy", "diff.add(ids[k][i])", "diff.add((k, ids[k][i]))"),
 (S, TS, "the donor's number always read from its row", "    if n is None:\n        n = contact_differences(rows[donor], human, columns)", "    if True:\n        n = contact_differences(rows[donor], human, columns)"),
 (S, TS, "the given number of differences ignored", "    if n is None:\n        n = contact_differences(rows[donor], human, columns)", "    n = 1 if n is not None else contact_differences(rows[donor], human, columns)"),
 (S, TS, "a gap of the table read from the alignment", "    if reason is not None:\n        return done(case_row(donor, complex_, gene, reason))", "    if False:\n        return done(case_row(donor, complex_, gene, reason))"),
 (S, TS, "the protein of the table not used", "    if ortholog is None:\n        donor_sub, placed, lost = u_subs(rows[donor])", "    if True:\n        donor_sub, placed, lost = u_subs(rows[donor])"),
 (S, TS, "no null said for a gene without an alignment", '    elif rows is None:\n        null = "not assessed: no alignment of the gene (the null draws from the alignments)"\n', ""),
 (S, TS, "the protein reading not marked", 'extra = {"reading": "protein", ', "extra = {"),
 (S, TS, "the alignment reading marked", "mt_substituted=len(mt_sub), unmodelled=unmodelled, pdb=pid, **extra)", 'mt_substituted=len(mt_sub), unmodelled=unmodelled, pdb=pid, reading="protein", **extra)'),
 (S, TS, "a gene without an alignment read from the protein as if it had columns", "ids = cols if cols is not None else {c: list(chains[c][\"num\"]) for c in u_idx}", "ids = {c: [None] * len(chains[c][\"seq\"]) for c in u_idx}"),
 (S, TS, "only the first copy given the donor's protein", "        for c in u_idx:\n            s_, p_, l_, placement[c] = protein_subs", "        for c in u_idx[:1]:\n            s_, p_, l_, placement[c] = protein_subs"),
]
only = sys.argv[1:]                       # names (or parts of names) of the breaks to run, all when none is given
breaks = [b for b in breaks if not only or any(o in b[2] for o in only)]
caught = 0
for name_file, tests, name, a, b in breaks:
    T = Path(name_file)
    orig = T.read_text()
    if orig.count(a) != 1:
        print("NOT APPLIED", name, flush=True); continue
    try:
        T.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "--no-sync", "pytest", "-q", "-x", "-p", "no:cacheprovider", *tests], capture_output=True, text=True)
    finally:
        T.write_text(orig)
    ok = r.returncode != 0
    caught += ok
    print("caught " if ok else "SURVIVED", name, flush=True)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) else 1)
