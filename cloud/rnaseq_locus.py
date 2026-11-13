#!/usr/bin/env python3
"""Summarise spliced RNA-seq alignments for the giant loci of one species (v2, lean).

Input on stdin: SAM records that samtools has ALREADY restricted to primary, mapped
alignments overlapping the locus spans (samtools view -F 0x904 -L loci.bed), so this
script only ever sees reads that touch a giant locus. Exon look-ups use bisect, coverage
uses difference arrays; the per-read cost is O(blocks * log exons).

Usage: ... | rnaseq_locus.py loci.json run_id total_primary_file > run.json
loci.json: [{"locus_id","md5","rep_accession","contig","strand","segments":[[s,e],...]}]
Per locus (MAPQ >= 10 counted as unique):
  reads_unique, reads_multi overlapping CDS exons or spliced inside the gene span
  introns: [start, end, exact_unique, exact_multimapped, alt_same_donor, alt_same_acceptor]
  exons:   [start, end, covered_frac, mean_depth] (unique reads)
  novel junctions (>= 2 unique reads) with an end inside a CDS exon, flag 1 if the
  other end leaves the gene span
"""
import json
import re
import sys
from bisect import bisect_right
from collections import Counter, defaultdict

loci = json.load(open(sys.argv[1]))
run_id = sys.argv[2]
total_file = sys.argv[3] if len(sys.argv) > 3 else ""
MINQ = 10
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

by_contig = defaultdict(list)
for L in loci:
    ex = sorted((int(a), int(b)) for a, b in L["segments"])
    L["_ex"] = ex
    L["_starts"] = [a for a, _ in ex]
    L["_lo"], L["_hi"] = ex[0][0], ex[-1][1]
    L["_introns"] = {(ex[k][1] + 1, ex[k + 1][0] - 1): k for k in range(len(ex) - 1)}
    L["_donor"] = {i[0]: k for i, k in L["_introns"].items()}
    L["_acc"] = {i[1]: k for i, k in L["_introns"].items()}
    L["_exact"], L["_exact_multi"] = Counter(), Counter()
    L["_alt_donor"], L["_alt_acc"], L["_novel"] = Counter(), Counter(), Counter()
    L["_diff"] = [[0] * (b - a + 2) for a, b in ex]
    L["_uniq"] = L["_multi"] = 0
    by_contig[L["contig"]].append(L)


def exon_index(L, pos):
    k = bisect_right(L["_starts"], pos) - 1
    return k if k >= 0 and pos <= L["_ex"][k][1] else None


n_in = 0
for line in sys.stdin:
    if line[0] == "@":
        continue
    f = line.split("\t", 6)
    flag = int(f[1])
    if flag & 0x904:
        continue
    Ls = by_contig.get(f[2])
    if not Ls:
        continue
    n_in += 1
    pos, mapq = int(f[3]), int(f[4])
    blocks, juncs, p = [], [], pos
    for n, op in CIG.findall(f[5]):
        n = int(n)
        if op in "M=X":
            blocks.append((p, p + n - 1))
            p += n
        elif op == "D":
            p += n
        elif op == "N":
            juncs.append((p, p + n - 1))
            p += n
    end = p - 1
    for L in Ls:
        if end < L["_lo"] - 2000 or pos > L["_hi"] + 2000:
            continue
        hits = []
        for a, b in blocks:                       # exons overlapped by each aligned block
            k = max(0, bisect_right(L["_starts"], a) - 1)
            while k < len(L["_ex"]) and L["_ex"][k][0] <= b:
                ea, eb = L["_ex"][k]
                if eb >= a:
                    hits.append((k, max(a, ea), min(b, eb)))
                k += 1
        if not hits and not juncs:
            continue
        if mapq < MINQ:
            L["_multi"] += 1
            for j in juncs:
                k = L["_introns"].get(j)
                if k is not None:
                    L["_exact_multi"][k] += 1
            continue
        L["_uniq"] += 1
        for k, lo, hi in hits:
            d = L["_diff"][k]
            ea = L["_ex"][k][0]
            d[lo - ea] += 1
            d[hi - ea + 1] -= 1
        for j in juncs:
            k = L["_introns"].get(j)
            if k is not None:
                L["_exact"][k] += 1
                continue
            if j[0] in L["_donor"]:
                L["_alt_donor"][L["_donor"][j[0]]] += 1
            if j[1] in L["_acc"]:
                L["_alt_acc"][L["_acc"][j[1]]] += 1
            if exon_index(L, j[0] - 1) is not None or exon_index(L, j[1] + 1) is not None:
                L["_novel"][j] += 1

total = None
if total_file:
    try:
        total = int(open(total_file).read().strip() or 0)
    except (OSError, ValueError):
        total = None
out = {"run": run_id, "primary_alignments": total, "primary_in_locus_spans": n_in,
       "primary_outside_loci": (total - n_in) if total is not None else None, "loci": []}
for L in loci:
    introns = [[a, b, L["_exact"][k], L["_exact_multi"][k], L["_alt_donor"][k], L["_alt_acc"][k]]
               for (a, b), k in sorted(L["_introns"].items())]
    exons = []
    for (a, b), d in zip(L["_ex"], L["_diff"]):
        n, run, covered, tot = b - a + 1, 0, 0, 0
        for x in range(n):
            run += d[x]
            if run:
                covered += 1
                tot += run
        exons.append([a, b, round(covered / n, 3), round(tot / n, 2)])
    novel = [[a, b, c, int(a - 1 < L["_lo"] or b + 1 > L["_hi"])]
             for (a, b), c in L["_novel"].most_common(50) if c >= 2]
    out["loci"].append({
        "locus_id": L["locus_id"], "md5": L["md5"], "rep_accession": L["rep_accession"],
        "reads_unique": L["_uniq"], "reads_multi": L["_multi"], "n_introns": len(introns),
        "introns_with_ge3_exact_reads": sum(1 for i in introns if i[2] >= 3),
        "introns_with_ge3_incl_multimapped": sum(1 for i in introns if i[2] + i[3] >= 3),
        "intron_columns": ["start", "end", "exact_unique", "exact_multimapped", "alt_same_donor",
                           "alt_same_acceptor"],
        "introns": introns, "exons": exons, "novel_junctions": novel,
        "novel_junctions_leaving_gene": sum(1 for x in novel if x[3]),
    })
json.dump(out, sys.stdout, separators=(",", ":"))
