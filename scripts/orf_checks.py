"""Cheap sequence-level checks on a single CDS region of a genomic contig.

Usage: orf_checks.py <fasta(.gz)> <contig_substring> <start> <end> <strand +|->

Reports: N count, GC, GC at codon positions 1/2/3 of the annotated frame, stop-codon
density per 1000 codons in the annotated frame and the two shifted frames on the
same strand and the three frames of the opposite strand, longest exact internal
repeat fraction (30-mers), and DNA tandem period.
A real protein-coding frame in a genome usually shows a distinct GC1/GC2/GC3
pattern (GC2 lowest; GC3 highest in GC-rich genomes) and many stops in the other
frames; a long ORF that exists only because all frames lack stops (GC-rich or
tandem-repeat DNA) does not.
"""
import gzip
import sys
from collections import Counter

fa, key, start, end, strand = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
op = gzip.open if fa.endswith(".gz") else open
seq, take = [], False
for line in op(fa, "rt"):
    if line.startswith(">"):
        take = key in line
        continue
    if take:
        seq.append(line.strip())
contig = "".join(seq).upper()
s = contig[start - 1:end]
comp = str.maketrans("ACGTN", "TGCAN")
if strand == "-":
    s = s.translate(comp)[::-1]
print(f"contig length {len(contig)}, region {start}-{end} ({len(s)} bp, strand {strand}), "
      f"distance to contig ends {start - 1} / {len(contig) - end}")
print(f"N: {s.count('N')}, GC: {(s.count('G') + s.count('C')) / len(s):.3f}, "
      f"contig GC: {(contig.count('G') + contig.count('C')) / len(contig):.3f}")
for pos in range(3):
    col = s[pos::3]
    print(f"  GC{pos + 1}: {(col.count('G') + col.count('C')) / len(col):.3f}")
STOPS = {"TAA", "TAG", "TGA"}


def stops(x, f):
    cod = [x[i:i + 3] for i in range(f, len(x) - 2, 3)]
    return 1000 * sum(c in STOPS for c in cod) / len(cod)


rc = s.translate(comp)[::-1]
print("stops per 1000 codons: annotated frame %.2f | same strand +1 %.2f, +2 %.2f | "
      "opposite strand %.2f, %.2f, %.2f" % (stops(s, 0), stops(s, 1), stops(s, 2),
                                            stops(rc, 0), stops(rc, 1), stops(rc, 2)))
k = 30
c = Counter(s[i:i + k] for i in range(len(s) - k + 1))
print(f"fraction of positions in repeated 30-mers: "
      f"{sum(1 for i in range(len(s) - k + 1) if c[s[i:i + k]] > 1) / len(s):.3f}")
best = sorted(((sum(1 for i in range(0, len(s) - p, 11) if s[i] == s[i + p]) /
                len(range(0, len(s) - p, 11)), p) for p in range(2, 600)), reverse=True)[:3]
print("DNA tandem periods (lag, identity):", [(p, round(f, 3)) for f, p in best])

# windowed frame signal: in 3,000-bp windows, is GC2 the lowest codon position?
W = 3000
ok = tot = 0
bad = []
for w in range(0, len(s) - W + 1, W):
    x = s[w:w + W]
    g = [(x[p::3].count("G") + x[p::3].count("C")) / len(x[p::3]) for p in range(3)]
    tot += 1
    if g[1] < g[0] and g[1] < g[2]:
        ok += 1
    else:
        bad.append(w // 3)
print(f"windows ({W} bp) with GC2 lowest: {ok}/{tot}; first failing windows at aa positions: {bad[:10]}")
