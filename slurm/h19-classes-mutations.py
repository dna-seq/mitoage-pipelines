"""Break each guard of H19's code (`mitoport classes`) in turn and check that a test fails.

Run by slurm/h19-classes.sbatch after the tests have passed on the unbroken code. Each mutation replaces one line of
the source, runs tests/test_classes.py (stopping at the first failure), and puts the file back byte for byte from a
copy (git is not installed on the compute nodes). Prints one line per mutation, "caught" or "SURVIVED", and exits
non-zero if any survived, was not applied, or the file did not come back.
"""
import hashlib
import signal
import subprocess
import sys
from pathlib import Path

CL = Path("src/mitoage_pipelines/mitoport/classes.py")
TESTS = ["tests/test_classes.py"]
MUTATIONS = [
    (CL, "the columns that touch a mitochondrial subunit are read",
     'if hrow[c] != "-" and c not in touching], dtype=int)', 'if hrow[c] != "-"], dtype=int)'),
    (CL, "the logarithm of time is not partialled out",
     "cov = [T[ok], np.log1p(T[ok])] + [stats.rankdata(e[ok]) for e in extra]",
     "cov = [T[ok]] + [stats.rankdata(e[ok]) for e in extra]"),
    (CL, "a covariate enters as its value, not its rank",
     "cov = [T[ok], np.log1p(T[ok])] + [stats.rankdata(e[ok]) for e in extra]",
     "cov = [T[ok], np.log1p(T[ok])] + [e[ok] for e in extra]"),
    (CL, "too few pairs are correlated anyway", "    if n < MIN_PAIRS:", "    if n < 2:"),
    (CL, "codons of different fourfold families share a site",
     "            fam = first == prefix", "            fam = np.isin(first, FOURFOLD_PREFIXES)"),
    (CL, "an unknown third base is a site",
     'purine, pyrimidine = np.isin(third, ["A", "G"]), np.isin(third, ["C", "T"])',
     'purine, pyrimidine = np.isin(third, ["A", "G", "N"]), np.isin(third, ["C", "T"])'),
    (CL, "two species with few shared fourfold sites get a distance",
     "return np.where(sites >= MIN_FOURFOLD, tv / np.maximum(sites, 1), np.nan)",
     "return np.where(sites >= 1, tv / np.maximum(sites, 1), np.nan)"),
    (CL, "a codon alignment out of frame is read", "        if A.shape[1] % 3:", "        if False:"),
    (CL, "a gene in several chains takes the chain farthest from the core",
     "            if best is None or cand[0] < best[0]:", "            if best is None or cand[0] > best[0]:"),
    (CL, "the mean distance, not the median",
     "cand = (float(np.median(d)), c, float(d.min()))", "cand = (float(np.mean(d)), c, float(d.min()))"),
    (CL, "only the reference chains are the core",
     "core = cKDTree(np.vstack([xyz[c] for c in sorted(xyz) if gene_of[c] in cx.mt]))",
     "core = cKDTree(np.vstack([xyz[c] for c in sorted(xyz) if c in cx.mt.values()]))"),
    (CL, "a species without a codon alignment is kept",
     "have = set.intersection(*(set(a) for a in mt_aln.values()), *(set(a) for a in codon_aln.values())) & set(tips)",
     "have = set.intersection(*(set(a) for a in mt_aln.values())) & set(tips)"),
    (CL, "a species that is no tip of the tree is kept",
     "have = set.intersection(*(set(a) for a in mt_aln.values()), *(set(a) for a in codon_aln.values())) & set(tips)",
     "have = set.intersection(*(set(a) for a in mt_aln.values()), *(set(a) for a in codon_aln.values()))"),
    (CL, "species of every class are read",
     'return sorted(set(species.slug[species["class"].eq(cls)]) & have & nuclear)',
     "return sorted(set(species.slug) & have & nuclear)"),
    (CL, "a protein with too few species is scored", "        if len(mine) < MIN_SPECIES:", "        if len(mine) < 2:"),
    (CL, "a class with too few species is scored", "    if len(sp) < MIN_SPECIES:", "    if len(sp) < 2:"),
    (CL, "the contact residues stand for the mitochondrial side",
     'sites[sites.gene.eq(g) & sites.group.eq("control")].aln_col.to_numpy() - 1 for g in cx.mt}',
     'sites[sites.gene.eq(g) & sites.group.eq("contact")].aln_col.to_numpy() - 1 for g in cx.mt}'),
    (CL, "the control columns are read one to the right",
     'sites[sites.gene.eq(g) & sites.group.eq("control")].aln_col.to_numpy() - 1 for g in cx.mt}',
     'sites[sites.gene.eq(g) & sites.group.eq("control")].aln_col.to_numpy() for g in cx.mt}'),
    (CL, "a gene's distance is not scaled before pooling", "        Ds.append(D / np.nanmean(D))", "        Ds.append(D)"),
    (CL, "the complex's own genes are among the other complexes'",
     "[g for g in all_cols if g not in control_cols[c]]) for c in D}", "list(all_cols)) for c in D}"),
    (CL, "the other complexes' genes are not partialled out",
     'row["given_other"] = pcorr(dn, dm, T, (do,))[0]', 'row["given_other"] = pcorr(dn, dm, T)[0]'),
    (CL, "the neutral distance is not partialled out of the protein coupling",
     'row["aa_given_tv4"] = pcorr(dn, dm, T, (s4,))[0]', 'row["aa_given_tv4"] = pcorr(dn, dm, T)[0]'),
    (CL, "the proteins are not partialled out of the neutral coupling",
     'row["tv4_given_aa"] = pcorr(dn, s4, T, (dm,))[0]', 'row["tv4_given_aa"] = pcorr(dn, s4, T)[0]'),
    (CL, "ranks are pooled over the complexes",
     "        out[m] = (stats.rankdata(v[m]) - 0.5) / m.sum()\n    return out",
     "        pass\n    return (stats.rankdata(v) - 0.5) / len(v)"),
    (CL, "the permutation mixes the complexes",
     "        for gidx in groups:\n            p[gidx] = ry[rng.permutation(gidx)]\n        hits += abs(np.corrcoef(rx, p)",
     "        p = ry[rng.permutation(len(ry))]\n        hits += abs(np.corrcoef(rx, p)"),
    (CL, "a complex with two subunits enters the stratified correlation",
     "    keep = _usable(strata, ok, MIN_STRATUM)", "    keep = _usable(strata, ok, 1)"),
    (CL, "the divergence is not partialled out of G2",
     "got = stratified_corr(own.median_dist, own.coupling, st, cov, seed=SEED + 10 * k + (cov is not None))",
     "got = stratified_corr(own.median_dist, own.coupling, st, None, seed=SEED + 10 * k + (cov is not None))"),
    (CL, "subunits with one to four interface columns count as touching",
     "got = stratified_diff(two.coupling, two.iface_cols >= MIN_IFACE, two.complex.to_numpy(), seed=SEED + 10 * k + 2)",
     "got = stratified_diff(own.coupling, own.iface_cols > 0, own.complex.to_numpy(), seed=SEED + 10 * k + 2)"),
    (CL, "one-sided group test", 'alternative="two-sided").pvalue)', 'alternative="greater").pvalue)'),
    (CL, "too few proteins compared anyway",
     "    if len(a) < MIN_PROTEINS or len(b) < MIN_PROTEINS:", "    if len(a) < 1 or len(b) < 1:"),
    (CL, "Holm over every readout",
     'res[res.status.eq("tested") & res.readout.isin(PRIMARY)].groupby("class")', 'res[res.status.eq("tested")].groupby("class")'),
    (CL, "an alignment without the human protein is placed on the structure",
     "interface_columns(model, chains[~chains.gene.isin(headless)], ct[ct.partner_is_mt == 0], nucaln)",
     "interface_columns(model, chains, ct[ct.partner_is_mt == 0], nucaln)"),
    (CL, "a missing contact map is not refused", '            if not (w / f).is_file():', '            if False:'),
]


def main() -> int:
    # Slurm's time limit (and kill) send SIGTERM: raise, so that `finally` puts the source back
    signal.signal(signal.SIGTERM, lambda *_: sys.exit("terminated: sources put back"))
    originals = {f: f.read_bytes() for f in {m[0] for m in MUTATIONS}}
    digests = {f: hashlib.sha256(b).hexdigest() for f, b in originals.items()}
    bad = 0
    try:
        for path, name, old, new in MUTATIONS:
            text = originals[path].decode()
            if text.count(old) != 1:
                print(f"NOT APPLIED (line found {text.count(old)} times): {name}"); bad += 1; continue
            path.write_text(text.replace(old, new, 1))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "--show-capture=no", *TESTS],
                               capture_output=True, text=True)
            failed = [l.split(" - ")[0].removeprefix("FAILED ") for l in r.stdout.splitlines() if l.startswith(("FAILED", "ERROR"))]
            if r.returncode == 0:
                print(f"SURVIVED: {name}"); bad += 1
            else:
                print(f"caught: {name}: {', '.join(failed) or r.stdout.strip().splitlines()[-1]}", flush=True)
            path.write_bytes(originals[path])
    finally:
        for f, b in originals.items():
            f.write_bytes(b)
    back = all(hashlib.sha256(f.read_bytes()).hexdigest() == d for f, d in digests.items())
    print(f"{len(MUTATIONS) - bad} of {len(MUTATIONS)} mutations caught; sources restored: {back}")
    return 0 if bad == 0 and back else 1


if __name__ == "__main__":
    sys.exit(main())
