#!/usr/bin/env python3
"""SAM (minimap2 --sam-hit-only) on stdin -> FASTQ of primary records on stdout.
Arguments: mate number appended to read names (1/2); optional cap on the number of reads
written. When the cap is reached the script exits, which closes the pipe and stops the
upstream download/decompression/alignment (the run's saturation rule).
Reverse-strand records are reverse-complemented back to the read orientation."""
import sys

mate = sys.argv[1] if len(sys.argv) > 1 else "1"
cap = int(sys.argv[2]) if len(sys.argv) > 2 else 0
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
w = sys.stdout.write
n = 0
for line in sys.stdin:
    if line[0] == "@":
        continue
    f = line.split("\t", 11)
    flag = int(f[1])
    if flag & 0x904:
        continue
    seq, qual = f[9], f[10]
    if flag & 0x10:
        seq, qual = seq.translate(COMP)[::-1], qual[::-1]
    if qual == "*":
        qual = "I" * len(seq)
    w(f"@{f[0]}/{mate}\n{seq}\n+\n{qual}\n")
    n += 1
    if cap and n >= cap:
        sys.stdout.flush()
        sys.stderr.write(f"sam2fq: cap of {cap} candidate reads reached\n")
        break
