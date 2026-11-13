"""Job 1, step 2: sequence-level checks on the genomic locus of each resolved giant
(from data/loci/specs.jsonl + data/loci/dna/*.fa.gz written by 10_locus_fetch.py).

Per sequence (all values measured):
  structure : n_exons, CDS length, mean/median/max exon, mean/max intron, genomic span
  assembly  : contig length, distance of the gene to both contig ends, N bases in the
              CDS, assembly-gap runs (>= 10 N) inside the gene span
  model     : translation of the exon chain vs the database protein (identity)
  frames    : internal stops in the annotated frame; stops per 1,000 codons in the
              other two frames on the same strand and the three opposite-strand frames
              (an ORF that is open in every frame is frame-independent)
  coding    : GC1/GC2/GC3; share of 3-kb CDS windows where GC2 is the lowest position
  repeats   : best DNA tandem lag (2..600 bp) and self-identity at that lag on the CDS
Output: data/locus_dna_checks.tsv
"""
import gzip
import json
import os
import re
import sys
from multiprocessing import Pool

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
LD = os.path.join(D, "loci")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG", "TGA"}
CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a + b + c] = aas[i]
            i += 1


def read_fasta_dict(path):
    out, name, buf = {}, None, []
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


# Fickett TESTCODE (Fickett 1982), tables as used in CPAT; >= 0.95 "coding", <= 0.74 "non-coding"
F_POS_PROB = {"A": [0.94, 0.68, 0.84, 0.93, 0.58, 0.68, 0.45, 0.34, 0.20, 0.22],
              "C": [0.80, 0.70, 0.70, 0.81, 0.66, 0.48, 0.51, 0.33, 0.30, 0.23],
              "G": [0.90, 0.88, 0.74, 0.64, 0.53, 0.48, 0.27, 0.16, 0.08, 0.08],
              "T": [0.97, 0.97, 0.91, 0.68, 0.69, 0.44, 0.54, 0.20, 0.09, 0.09]}
F_POS_W = {"A": 0.26, "C": 0.18, "G": 0.31, "T": 0.33}
F_POS_PARA = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
F_CON_PROB = {"A": [0.28, 0.49, 0.44, 0.55, 0.62, 0.49, 0.67, 0.65, 0.81, 0.21],
              "C": [0.82, 0.64, 0.51, 0.64, 0.59, 0.59, 0.43, 0.44, 0.39, 0.31],
              "G": [0.40, 0.54, 0.47, 0.64, 0.64, 0.73, 0.41, 0.41, 0.33, 0.29],
              "T": [0.28, 0.24, 0.39, 0.40, 0.55, 0.75, 0.56, 0.69, 0.51, 0.58]}
F_CON_W = {"A": 0.11, "C": 0.12, "G": 0.15, "T": 0.14}
F_CON_PARA = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.17, 0]


def fickett(x):
    n = len(x)
    score = 0.0
    for b in "ACGT":
        c = [x[p::3].count(b) for p in range(3)]
        pos = max(c) / (min(c) + 1.0)
        for i, v in enumerate(F_POS_PARA):
            if pos >= v:
                score += F_POS_PROB[b][i] * F_POS_W[b]
                break
        con = x.count(b) / n
        for i, v in enumerate(F_CON_PARA):
            if con >= v:
                score += F_CON_PROB[b][i] * F_CON_W[b]
                break
    return score


def stops_per_k(x, f):
    cod = [x[i:i + 3] for i in range(f, len(x) - 2, 3)]
    return 1000 * sum(c in STOPS for c in cod) / max(1, len(cod))


def best_lag(s, maxlag=600, step=11):
    best = (0.0, 0)
    for p in range(2, min(maxlag, len(s) // 2)):
        idx = range(0, len(s) - p, step)
        m = sum(1 for i in idx if s[i] == s[i + p]) / max(1, len(idx))
        if m > best[0]:
            best = (m, p)
    return best


def analyze(args):
    spec, region, protein = args
    out = {"md5": spec["md5"], "rep_accession": spec["rep_accession"], "locus_source": spec["source"],
           "contig": spec["contig"], "strand": spec["strand"], "n_contigs_in_model": spec["n_contigs_in_model"],
           "locus_note": spec.get("note", "")}
    if not spec["segments"] or region is None:
        out["locus_note"] = (out["locus_note"] + "; no DNA").strip("; ")
        return out
    rstart, dna, clen = region
    segs = spec["segments"]
    seqs = []
    for a, b in segs:
        piece = dna[a - rstart:b - rstart + 1].upper()
        seqs.append(piece.translate(COMP)[::-1] if spec["strand"] == "-" else piece)
    cds = "".join(seqs)
    lo, hi = min(a for a, b in segs), max(b for a, b in segs)
    span = dna[lo - rstart:hi - rstart + 1].upper()
    exl = [b - a + 1 for a, b in segs]
    ordered = sorted(segs)
    introns = [ordered[k + 1][0] - ordered[k][1] - 1 for k in range(len(ordered) - 1)]
    exs = sorted(exl)
    out.update({
        "n_exons": len(segs), "cds_bp": len(cds), "mean_exon_bp": round(sum(exl) / len(exl)),
        "median_exon_bp": exs[len(exs) // 2], "max_exon_bp": max(exl),
        "mean_intron_bp": round(sum(introns) / len(introns)) if introns else 0,
        "max_intron_bp": max(introns) if introns else 0, "genomic_span_bp": hi - lo + 1,
        "contig_length": clen if clen else "",
        "dist_to_contig_start": lo - 1, "dist_to_contig_end": (clen - hi) if clen else "",
        "N_in_cds": cds.count("N"), "gap_runs_in_span": len(re.findall(r"N{10,}", span)),
        "N_in_span": span.count("N"),
        "gc_cds": round((cds.count("G") + cds.count("C")) / max(1, len(cds)), 3),
    })
    # translation of the exon chain vs the database protein
    tr = "".join(CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3)).rstrip("*")
    n = min(len(tr), len(protein))
    same = sum(1 for i in range(1, n) if tr[i] == protein[i])  # skip codon 1 (alt starts)
    out["translation_identity"] = round(same / max(1, n - 1), 4)
    out["translation_len_diff"] = len(tr) - len(protein)
    for p in range(3):
        col = cds[p::3]
        out[f"gc{p + 1}"] = round((col.count("G") + col.count("C")) / max(1, len(col)), 3)
    W, ok, tot = 3000, 0, 0
    for w in range(0, len(cds) - W + 1, W):
        x = cds[w:w + W]
        g = [(x[q::3].count("G") + x[q::3].count("C")) for q in range(3)]
        tot += 1
        ok += g[1] < g[0] and g[1] < g[2]
    out["coding_windows"] = tot
    out["coding_windows_gc2_lowest_frac"] = round(ok / tot, 3) if tot else ""
    fk = [fickett(cds[w:w + 300]) for w in range(0, len(cds) - 299, 300)]
    out["fickett_windows_300bp"] = len(fk)
    out["fickett_coding_frac"] = round(sum(1 for v in fk if v >= 0.95) / len(fk), 3) if fk else ""
    out["fickett_noncoding_frac"] = round(sum(1 for v in fk if v <= 0.74) / len(fk), 3) if fk else ""
    out["fickett_median"] = round(sorted(fk)[len(fk) // 2], 3) if fk else ""
    rc = cds.translate(COMP)[::-1]
    other = [stops_per_k(cds, 1), stops_per_k(cds, 2), stops_per_k(rc, 0), stops_per_k(rc, 1), stops_per_k(rc, 2)]
    out["internal_stops_annot_frame"] = sum(1 for i in range(0, len(cds) - 5, 3) if cds[i:i + 3] in STOPS)
    out["min_stops_per_kcodon_other_frames"] = round(min(other), 2)
    out["n_other_frames_open"] = sum(1 for v in other if v < 5)   # of 5; < 5 stops per 1,000 codons
    out["mean_stops_per_kcodon_other_frames"] = round(sum(other) / 5, 2)
    m, p = best_lag(cds)
    out["dna_tandem_lag_bp"], out["dna_tandem_identity"] = p, round(m, 3)
    return out


if __name__ == "__main__":
    specs = [json.loads(l) for l in open(os.path.join(LD, "specs.jsonl"))]
    extra = os.path.join(LD, "control_specs.jsonl")
    if os.path.exists(extra):
        specs += [json.loads(l) for l in open(extra)]
    index = json.load(open(os.path.join(LD, "dna_index.json")))
    clen = json.load(open(os.path.join(LD, "contig_lengths.json")))
    prot = read_fasta_dict(os.path.join(D, "fasta", "unique_ge35992.fasta"))
    cpath = os.path.join(D, "fasta", "controls.fasta")
    if os.path.exists(cpath):
        prot.update(read_fasta_dict(cpath))
    by_md5 = {}
    for lid, v in index.items():
        for m in v["md5s"]:
            by_md5[m] = lid
    jobs = []
    cache = {}
    for s in specs:
        lid = by_md5.get(s["md5"])
        region = None
        if lid and index[lid]["how"] not in ("fail", "fail-too-large"):
            if lid not in cache:
                fa = read_fasta_dict(os.path.join(LD, "dna", lid + ".fa.gz"))
                cache.clear()
                cache[lid] = next(iter(fa.values()))
            region = (index[lid]["start"], cache[lid], clen.get(s["contig"]))
        jobs.append((s, region, prot.get(s["md5"], "")))
    with Pool(int(os.environ.get("NPROC", "8"))) as pool:
        res = pool.map(analyze, jobs, chunksize=1)
    cols = []
    for r in res:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(os.path.join(D, "locus_dna_checks.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in res:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print("rows", len(res), "with DNA", sum(1 for r in res if "cds_bp" in r))
