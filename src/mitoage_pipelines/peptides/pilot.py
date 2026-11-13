"""Pilot: are the published mtDNA-encoded peptides' reading frames kept across species?

For each peptide: find it in the human row of the gene's alignment (this also checks the
sequence typed here against the data), project its columns to every species, and test the frame.
Null: every other window of the same length in the same gene that is stop-free in human.

Writes one table per peptide and summary.txt.
"""
import random
import statistics
from pathlib import Path

import numpy as np
import polars as pl
from Bio.Data import CodonTable

from ..seqio import read_fasta


def table(n, agr_arg=False):
    t = CodonTable.unambiguous_dna_by_id[n]
    fwd = dict(t.forward_table)
    stops = set(t.stop_codons)
    if agr_arg:
        for c in ("AGA", "AGG"):
            stops.discard(c); fwd[c] = "R"
    return fwd, stops


def translate(nt, code):
    fwd, stops = code
    out = []
    for i in range(0, len(nt) - 2, 3):
        c = nt[i:i + 3]
        if c in stops: out.append("*")
        else: out.append(fwd.get(c, "X"))
    return "".join(out)


STD = table(1)

# name, file, peptide (a number = search by length), code used for human / frames to search
PEPTIDES = [
    ("humanin", "rRNA-16S.nt", "MAPRGFSCLLLLTSEIDLPVKRRA", "std"),
    ("SHLP6", "rRNA-16S.nt", "MLDQDIPMVQPLLKVRLFND", "std"),
    ("MOTS-c", "rRNA-12S.nt", "MRWQEMGYIFYPRKLR", "std"),
    ("MTALTND4", "ND4.codon", 99, "any"),
    ("MTALTCO1", "COX1.codon", 259, "any"),
]


def human_orfs(h, want, mode):
    """(start, codons, code name) of ORFs in the ungapped human row matching `want`."""
    found = []
    codes = {"std": STD, "mito2": table(2), "mito2_AGR=R": table(2, True)}
    if mode == "std": codes = {"std": STD}
    for cname, code in codes.items():
        for frame in range(3):
            prot = translate(h[frame:], code)
            if isinstance(want, str):
                i = prot.find(want)
                while i >= 0:
                    found.append((frame + 3 * i, len(want), cname)); i = prot.find(want, i + 1)
            else:
                pos = 0
                for seg in prot.split("*"):
                    m = seg.find("M")
                    while m >= 0:
                        if len(seg) - m == want:
                            found.append((frame + 3 * (pos + m), want, cname))
                        m = seg.find("M", m + 1)
                    pos += len(seg) + 1
    return found


def columns(row):
    """alignment column of each ungapped position"""
    return [i for i, c in enumerate(row) if c != "-"]


def window_state(seg, n, code):
    """seg: species' ungapped nt in the window. -> (frame kept, codons before first stop)"""
    prot = translate(seg, code)
    stop = prot.find("*")
    open_len = len(prot) if stop < 0 else stop
    return (len(seg) - 3 * n) % 3 == 0, open_len


def main(data: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    species = pl.read_parquet(data / "species.parquet").select(
        "slug", "class", "order", "max_longevity_years", "transl_table")
    META = {r["slug"]: r for r in species.iter_rows(named=True)}

    report = []
    for name, fname, want, mode in PEPTIDES:
        aln = {k: v.upper() for k, v in read_fasta(data / "alignments" / f"{fname}.fasta.gz").items()}
        hrow = aln["homo-sapiens"]; h = hrow.replace("-", ""); hcols = columns(hrow)
        hits = human_orfs(h, want, mode)
        print(f"## {name}: human hits {[(s, n, c, translate(h[s:s+3*n+3], table(2, True) if 'AGR' in c else (STD if c=='std' else table(2)))) for s, n, c in hits]}")
        if not hits:
            report.append(f"{name}\tNOT FOUND in the human row"); continue
        # prefer a hit in a frame other than the host's for the codon files
        hits.sort(key=lambda x: (x[0] % 3 == 0 and fname.endswith("codon"), x[2]))
        start, n, cname = hits[0]
        human_code = {"std": STD, "mito2": table(2), "mito2_AGR=R": table(2, True)}[cname]

        def scan(s):
            """per species state for the window starting at human position s"""
            c0, c1 = hcols[s], hcols[s + 3 * n - 1] + 1
            res = {}
            for slug, row in aln.items():
                seg = row[c0:c1].replace("-", "")
                if len(seg) < 3 * n * 0.5: res[slug] = None; continue   # mostly gap here
                m = META.get(slug)
                if cname == "std": code = STD
                else: code = table(m["transl_table"] if m else 2, "AGR" in cname)
                kept_frame, open_len = window_state(seg, n, code)
                res[slug] = (seg[:3], kept_frame, open_len, len(seg) // 3, translate(seg, code))
            return res

        def open_fraction(res, slugs):
            ok = tot = 0
            for s in slugs:
                r = res.get(s)
                if r is None: continue
                tot += 1; ok += r[2] >= 0.9 * min(n, r[3])
            return ok / tot if tot else float("nan"), tot

        res = scan(start)
        classes = {}
        for slug in aln:
            m = META.get(slug)
            if m: classes.setdefault(m["class"], []).append(slug)
        mammals = classes.get("Mammalia", [])

        # null: other windows of the same length, same frame logic, stop-free in human
        rng = random.Random(1)
        cand = [s for s in range(0, len(h) - 3 * n - 3) if abs(s - start) >= 3 * n or (s - start) % 3]
        cand = [s for s in cand if "*" not in translate(h[s:s + 3 * n], human_code)]
        if fname.endswith("codon"): cand = [s for s in cand if s % 3 == start % 3]  # same alternative frame
        rng.shuffle(cand); cand = cand[:300]
        null = [open_fraction(scan(s), mammals)[0] for s in cand]
        obs, ntot = open_fraction(res, mammals)
        pct = sum(x >= obs for x in null) / len(null) if null else float("nan")

        lines = [f"\n# {name} ({fname}; human start {start}, {n} codons, read with {cname})"]
        lines.append(f"mammals: open in {obs:.3f} of {ntot}; null windows (stop-free in human, n={len(null)}): "
                     f"median {statistics.median(null) if null else float('nan'):.3f}, share of null >= observed {pct:.3f}")
        lines.append("class\tspecies\twindow_present\tframe_kept\topen>=90%\tATG_start\topen_and_ATG\tmedian_identity_to_human")
        hpep = translate(h[start:start + 3 * n], human_code)
        with open(out / f"{name}.tsv", "w") as f:
            f.write("slug\tclass\torder\tmax_longevity_years\tstart_codon\tframe_kept\topen_codons\twindow_codons\tidentity\tpeptide\n")
            for cls, slugs in sorted(classes.items(), key=lambda kv: -len(kv[1])):
                pres = fk = op = atg = both = 0; ids = []
                for s in slugs:
                    r = res.get(s)
                    if r is None: continue
                    pres += 1; fk += r[1]
                    o = r[2] >= 0.9 * min(n, r[3]); a = r[0] == "ATG"
                    op += o; atg += a; both += o and a
                    ident = sum(x == y for x, y in zip(r[4], hpep)) / n if r[1] and r[3] == n else float("nan")
                    if ident == ident: ids.append(ident)
                    m = META[s]
                    f.write(f"{s}\t{cls}\t{m['order']}\t{m['max_longevity_years']}\t{r[0]}\t{int(r[1])}\t{r[2]}\t{r[3]}\t{ident:.3f}\t{r[4]}\n")
                if len(slugs) >= 20:
                    lines.append(f"{cls}\t{len(slugs)}\t{pres}\t{fk}\t{op}\t{atg}\t{both}\t{statistics.median(ids) if ids else float('nan'):.2f}")
        # lifespan, mammals only, not corrected for the tree
        kept = [np.log10(META[s]["max_longevity_years"]) for s in mammals if res.get(s) and res[s][2] >= 0.9 * min(n, res[s][3]) and res[s][0] == "ATG" and META[s]["max_longevity_years"]]
        lost = [np.log10(META[s]["max_longevity_years"]) for s in mammals if res.get(s) and not (res[s][2] >= 0.9 * min(n, res[s][3]) and res[s][0] == "ATG") and META[s]["max_longevity_years"]]
        if kept and lost:
            lines.append(f"mammals, raw (no tree): kept n={len(kept)} median lifespan {10**statistics.median(kept):.1f} y; lost n={len(lost)} median {10**statistics.median(lost):.1f} y")
        report += lines
        print("\n".join(lines), flush=True)

    (out / "summary.txt").write_text("\n".join(report) + "\n")
