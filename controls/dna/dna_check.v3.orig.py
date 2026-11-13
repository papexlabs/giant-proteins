#!/usr/bin/env python3
"""Read-level check of a giant single-ORF gene in its own assembly (DNA re-mapping).

Job (JSON file): {"job_id", "assembly_url" | "assembly_path", "contig", "gene_start", "gene_end",
  "strand", "illumina": [url_or_path_1, url_or_path_2] (paired), "ont": [url_or_path, ...]}

Steps
  1. assembly -> reference (all contigs: the rest of the assembly is the decoy)
  2. stream Illumina pairs (curl, named pipes, no storage) into minimap2 -ax sr against the
     reference, keep mapped reads only (--sam-hit-only) -> BAM; ONT likewise with map-ont
  3. negative control: the same reads that touch the target contig are re-mapped to a copy of
     the assembly where the target contig is cut in the middle of the gene and its second
     half is reverse-complemented (an in-silico chimeric join that does not exist in the
     genome); the same metrics must FAIL at that join
  4. metrics over the gene (1-kb windows), using primary alignments with MAPQ >= 20 and
     identity >= 97 % (NM / aligned length <= 0.03, which removes cross-mapping reads from
     related organisms in a metagenome):
       depth per window; share of gene bases at depth >= 1 and >= 5
       physical coverage: proper pairs (same contig, FR, insert <= 3x median) whose fragment
         spans each window boundary by >= 20 bp on both sides; ONT reads spanning each
         boundary by >= 500 bp on both sides
       clip pile-ups: positions where >= max(5, 20 % of depth) reads are soft-clipped >= 20 bp
       discordant pairs per window (mate on another contig, wrong orientation, insert > 5x median)
       gene median depth / contig median depth outside the gene
       ORF integrity on the reads' consensus: CDS positions (depth >= 5) where >= 50 % of
         reads carry an indel, or a substitution that creates a stop codon in frame
  5. verdict (thresholds in VERDICT below) -> <out_dir>/<job_id>.json

Rules v3 (2026-11-12). Up to v2 the only repeat mask was read depth per 1-kb window (unique
depth < 10 % of all-MAPQ depth), and every other boundary counted as single-copy. That misses a
repeat shorter than a window or one that straddles a window edge: on two Omnitrophota MAGs a
1.23-kb internal repeat (two copies, 99.8 % identical, 10.4 kb apart) sat across a window edge,
half of each window stayed unique, and the boundary inside the repeat got "0 unique spanning
pairs -> against", a count that no fragment of that library could ever exceed. Since v3 the
break test asks the library first: error-free pairs with the library's own insert and read
length are simulated from the assembly every 2 bp and mapped back the same way as the reads;
u(x) is the share placed uniquely (MAPQ >= 20, true place) across x, and E(x) = u(x) x the
gene's median spanning pairs is what the reads should show at x if the assembly is right.
  - a boundary is tested (break if 0 pairs; "weakest boundary" figure) only where E >= 5, the
    same minimum already required of an observed boundary
  - stretches with E < 5 are blocks this library cannot place; each must be CROSSED as a whole
    by >= 5 unique pairs or >= 3 long reads (500-bp flanks); otherwise it is an unresolved
    repeat and the best result is "supported outside N kb of repeat" (grade unchanged)
No threshold changed. Short repeats between window edges, ignored before, now need crossing.

Usage: dna_check.py job.json out_dir [--threads N] [--max-pairs N] [--assembly local.fa.gz]
       dna_check.py job.json out_dir --reverdict old_result.json [--assembly local.fa.gz]
         [--read-len N] (re-applies the current rules to a stored result without reads: the
          mappability is simulated again; crossing of a block is then known only when no accepted
          fragment can be that long or a boundary inside it has < 5 pairs, else it counts as not
          crossed)
"""
import argparse
from bisect import bisect_left, bisect_right
import gzip
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

UA = "PapexResearchBot/1.0 (+https://papex.org)"
ap = argparse.ArgumentParser()
ap.add_argument("job")
ap.add_argument("out")
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--max-pairs", type=int, default=0, help="cap read pairs (testing / subsample)")
ap.add_argument("--work", default="")
ap.add_argument("--assembly", default="", help="local assembly FASTA (instead of assembly_url)")
ap.add_argument("--reverdict", default="", help="stored result JSON: re-apply the current rules, no reads")
ap.add_argument("--read-len", type=int, default=0, help="Illumina read length for --reverdict (else 150)")
args = ap.parse_args()
RULES = "v3 2026-11-12: break test only where the library can place pairs (simulated from the assembly)"
job = json.load(open(args.job))
JID = job["job_id"]
W = args.work or (tempfile.mkdtemp(prefix="dnacheck_") if args.reverdict else os.path.join(args.out, "tmp_" + JID))
os.makedirs(W, exist_ok=True)
os.makedirs(args.out, exist_ok=True)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG", "TGA"}
MINQ, MAX_DIV, WIN = 20, 0.03, 1000          # MAX_DIV for Illumina; long reads: 0.15
VERDICT = {"min_frac_depth5": 0.99, "min_pairs_per_boundary": 5, "min_ont_per_boundary": 3,
           "max_discordant_frac": 0.10, "depth_ratio": (0.5, 2.0)}


def log(*a):
    print(time.strftime("%H:%M:%S"), JID, *a, flush=True)


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def src(x):
    """curl stream for URLs, cat for local files (gzip handled by minimap2)."""
    if x.startswith(("http://", "https://")):
        return f"curl -sfL --connect-timeout 30 -A '{UA}' '{x}'"
    return f"cat '{x}'"


def read_fasta(path):
    seqs, name, buf = {}, None, []
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


SIM_STEP, SIM_BIN = 2, 10


def library_mappability(ref_path, ctg, ctg_seq, lo, hi, insert, read_len):
    """Where this library could show a break at all if the assembly were right. Error-free FR
    pairs with the library's median insert and read length are cut from the assembly every 2 bp
    across ctg:lo-hi and mapped back exactly as the real reads were (minimap2 -ax sr against the
    whole assembly). u(x) = share of the simulated fragments covering x by >= 20 bp on both
    sides whose left read is a proper pair with MAPQ >= 20 at its true place. Returns u per
    10-bp bin (bin i starts at lo + 10 i); None where no fragment can be cut (Ns)."""
    insert = int(round(insert))
    rl = min(int(read_len), insert)
    n = hi - lo + 1
    tot, good = [0] * (n + 1), [0] * (n + 1)

    def add(arr, x, y):
        x, y = max(x, lo), min(y, hi)
        if x <= y:
            arr[x - lo] += 1
            arr[y - lo + 1] -= 1
    r1p, r2p = os.path.join(W, "sim_1.fq"), os.path.join(W, "sim_2.fq")
    qual = "I" * rl
    with open(r1p, "w") as f1, open(r2p, "w") as f2:
        for st in range(max(1, lo - insert + 21), min(len(ctg_seq) - insert + 1, hi - 20) + 1, SIM_STEP):
            frag = ctg_seq[st - 1:st - 1 + insert]
            if "N" in frag:
                continue
            f1.write(f"@f{st}\n{frag[:rl]}\n+\n{qual}\n")
            f2.write(f"@f{st}\n{frag[-rl:].translate(COMP)[::-1]}\n+\n{qual}\n")
            add(tot, st + 20, st + insert - 21)
    p = subprocess.Popen(f"minimap2 -t {args.threads} -ax sr {ref_path} {r1p} {r2p} 2>/dev/null "
                         f"| samtools view -F 0x904 -f 0x40 -", shell=True, stdout=subprocess.PIPE, text=True,
                         executable="/bin/bash")
    for line in p.stdout:
        f = line.split("\t", 9)
        flag, st = int(f[1]), int(f[0][1:])
        if flag & 0x2 and not flag & 0x10 and int(f[4]) >= MINQ and f[2] == ctg and int(f[3]) == st \
                and int(f[8]) == insert:
            add(good, st + 20, st + insert - 21)
    p.wait()
    u, ct, cg = [], 0, 0
    for i in range(n):
        ct += tot[i]
        cg += good[i]
        if i % SIM_BIN == 0:
            u.append(round(cg / ct, 3) if ct else None)
    return u


def spanning_block(frags, s, e, m=500):
    return sum(1 for a, b in frags if a <= s - m and b >= e + m)


def union_kb(ivs):
    tot, cur = 0, None
    for s, e in sorted(ivs):
        if cur and s <= cur[1] + 1:
            cur[1] = max(cur[1], e)
        else:
            if cur:
                tot += cur[1] - cur[0] + 1
            cur = [s, e]
    if cur:
        tot += cur[1] - cur[0] + 1
    return round(tot / 1000, 2)


def classify(res, kind, lo, u=None, frags=None, sim=None):
    """Boundary classes and the break test, from the stored window arrays (so a stored result
    can be re-judged). Illumina: E(x) = u(x) x S is the number of unique spanning pairs the
    library is expected to give at x if the assembly is right (S = the gene's median spanning
    pairs per boundary). A boundary is testable only if E >= 5 (the same minimum as for an
    observed boundary) and neither window is a read-depth repeat window; the break test and the
    weakest-boundary figure use testable boundaries only. Stretches with E < 5 are blocks this
    library cannot place; each must be crossed as a whole (>= 5 unique pairs >= 20 bp beyond each
    end, or >= 3 long reads >= 500 bp beyond each end) or it stays unresolved."""
    wq, wa, span = res["window_depth_unique"], res["window_depth_all_mapq"], res["spanning_per_boundary"]
    hi = lo + len(wq) * WIN - 1
    bounds = [lo + WIN * i for i in range(1, len(wq))]
    rep = [wa[i] >= 10 and wq[i] < 0.1 * wa[i] for i in range(len(wa))]
    need_depth = 10 if kind == "illumina" else 5
    margin = 20 if kind == "illumina" else 500
    need = VERDICT["min_pairs_per_boundary"] if kind == "illumina" else VERDICT["min_ont_per_boundary"]
    S = res.get("median_spanning_per_boundary") or 0
    low_bins = set()
    if kind == "illumina" and u is not None:
        low_bins = {i for i, v in enumerate(u) if (v or 0) * S < need}
    untestable = [((x - lo) // SIM_BIN) in low_bins for x in bounds]
    if kind == "illumina" and u is not None:
        res["expected_spanning_per_boundary"] = [round((u[(x - lo) // SIM_BIN] or 0) * S, 1)
                                                 if (x - lo) // SIM_BIN < len(u) else None for x in bounds]
    single = [k for k in range(len(bounds)) if not rep[k] and not rep[k + 1] and not untestable[k]]
    res["boundaries_untestable_for_library"] = [bounds[k] for k in range(len(bounds)) if untestable[k]]
    res["min_spanning_unique_boundaries"] = min((span[k] for k in single), default=None)
    res["break_boundaries"] = [bounds[k] for k in single if span[k] == 0 and min(wq[k], wq[k + 1]) >= need_depth]
    blocks, i, nb = [], 0, (len(u) if u else 0)
    while i < nb:
        if i in low_bins:
            j = i
            while j + 1 < nb and (j + 1) in low_bins:
                j += 1
            bs, be = lo + i * SIM_BIN, min(hi, lo + j * SIM_BIN + SIM_BIN - 1)
            vals = [u[k] for k in range(i, j + 1) if u[k] is not None]
            bb = {"start": bs, "end": be, "bp": be - bs + 1, "min_u": min(vals) if vals else None}
            inside = [k for k in range(len(bounds)) if bs <= bounds[k] <= be]
            bb["boundaries_inside"] = [bounds[k] for k in inside]
            bb["spanning_at_boundaries_inside"] = [span[k] for k in inside]
            longest = 3 * res["median_insert"] if res.get("median_insert") else None
            if frags is not None:
                bb["crossing"], bb["crossing_how"] = spanning_block(frags, bs, be, margin), "counted"
            elif longest and bb["bp"] + 2 * margin > longest:
                bb["crossing"] = 0
                bb["crossing_how"] = (f"none possible: block + 2 x {margin} bp > {longest:.0f} bp, the longest "
                                      "fragment the check accepts (3 x median insert)")
            elif inside and min(span[k] for k in inside) < need:
                bb["crossing"] = min(span[k] for k in inside)
                bb["crossing_how"] = "at most this many: fewest pairs at a boundary inside the block"
            else:
                bb["crossing"], bb["crossing_how"] = None, "not measured (needs the reads); counted as not crossed"
            bb["crossed"] = bb["crossing_how"] == "counted" and bb["crossing"] >= need
            blocks.append(bb)
            i = j + 1
        else:
            i += 1
    res["library_unplaceable_blocks"] = blocks[:50]
    res["library_unplaceable_blocks_n"] = len(blocks)
    res["unplaceable_kb_not_crossed"] = union_kb([(b["start"], b["end"]) for b in blocks if not b["crossed"]])
    if sim:
        sim = dict(sim)
        sim["bins"] = nb
        sim["frac_bins_testable"] = round(1 - len(low_bins) / nb, 4) if nb else None
        res["library_mappability"] = sim
    res["rules"] = RULES
    return res


# ---- 1. reference --------------------------------------------------------------------
t0 = time.time()
asm = os.path.join(W, "assembly.fa.gz")
if not os.path.exists(asm):
    if args.assembly:
        shutil.copy(args.assembly, asm)
    elif job.get("assembly_path"):
        shutil.copy(job["assembly_path"], asm)
    else:
        sh(f"curl -sfL --retry 5 -A '{UA}' -o {asm} '{job['assembly_url']}'")
seqs = read_fasta(asm)
# ENA headers look like ENA|DAOTDT010000001|DAOTDT010000001.1 ...: normalise to the accession
norm = {}
for k, v in seqs.items():
    acc = k.split("|")[-1] if "|" in k else k
    norm[acc] = v
seqs = norm
contig = job["contig"]
if contig not in seqs:
    alt = [k for k in seqs if k.split(".")[0] == contig.split(".")[0]]
    if not alt:
        sys.exit(f"contig {contig} not in assembly")
    contig = alt[0]
gs, ge = int(job["gene_start"]), int(job["gene_end"])
ref = os.path.join(W, "ref.fa")
with open(ref, "w") as f:
    for k, v in seqs.items():
        f.write(f">{k}\n{v}\n")
# negative control: cut the target contig at the gene midpoint, reverse-complement the 2nd part
mid = (gs + ge) // 2
chim = os.path.join(W, "ref_chimera.fa")
with open(chim, "w") as f:
    for k, v in seqs.items():
        if k == contig:
            v = v[:mid] + v[mid:].translate(COMP)[::-1]
        f.write(f">{k}\n{v}\n")
clen = len(seqs[contig])
log(f"assembly {len(seqs)} contigs, target {contig} {clen} bp, gene {gs}-{ge} ({ge-gs+1} bp)")
CTL = (mid - 5000, mid + 5000)        # the negative control looks at the junction neighbourhood
CHIM_SEQ = seqs[contig][:mid] + seqs[contig][mid:].translate(COMP)[::-1]


def judge_illumina(m, lo, hi, ref_path, ctg_seq, read_len, frags=None):
    """Library mappability for this region, then the boundary classes (classify)."""
    t = time.time()
    u = library_mappability(ref_path, contig, ctg_seq, lo, hi, m["median_insert"], read_len)
    sim = {"insert": int(round(m["median_insert"])), "read_length": int(read_len), "step_bp": SIM_STEP,
           "bin_bp": SIM_BIN, "seconds": round(time.time() - t)}
    return classify(m, "illumina", lo, u, frags, sim)

# ---- 2. mapping (streamed) ------------------------------------------------------------
cap = f"| head -n {4 * args.max_pairs}" if args.max_pairs else ""
def processed(logf):
    try:
        return sum(int(m) for m in re.findall(r"mapped (\d+) sequences", open(logf).read()))
    except OSError:
        return 0


def map_with_check(kind, cmd, logf, out_bam, expected):
    """Run a streamed mapping; a dropped download ends the stream early, so the number of
    reads minimap2 saw is compared with the run's read count (ENA) and the whole mapping is
    redone (from byte 0, into a fresh BAM: nothing is counted twice) up to 3 attempts."""
    for attempt in range(1, 4):
        if os.path.exists(out_bam + ".ok"):
            return
        sh(cmd)
        n = processed(logf)
        if args.max_pairs or not expected or n >= 0.995 * expected:
            sh(f"samtools index {out_bam} && touch {out_bam}.ok")
            log(f"{kind} mapped {n} reads (expected {expected}) {time.time()-t0:.0f}s")
            return
        log(f"{kind} stream incomplete: {n} of {expected} reads, attempt {attempt}; redoing")
    sys.exit(f"{kind}: stream kept failing")


bam = os.path.join(W, "illumina.bam")
if job.get("illumina") and not args.reverdict:
    r1, r2 = job["illumina"]
    map_with_check("illumina", f"minimap2 -t {args.threads} -ax sr --sam-hit-only {ref} "
                   f"<({src(r1)} | gzip -dc {cap}) <({src(r2)} | gzip -dc {cap}) 2> {W}/mm_illumina.log "
                   f"| samtools sort -@4 -o {bam} -", f"{W}/mm_illumina.log", bam, job.get("illumina_expected_reads"))
ont_bam = os.path.join(W, "ont.bam")
if job.get("ont") and not args.reverdict:
    files = " ".join(f"<({src(u)} | gzip -dc)" for u in job["ont"])
    map_with_check("ont", f"minimap2 -t {args.threads} -ax map-ont --sam-hit-only {ref} {files} 2> {W}/mm_ont.log "
                   f"| samtools sort -@4 -o {ont_bam} -", f"{W}/mm_ont.log", ont_bam, job.get("ont_expected_reads"))


def reads_on_contig_to_fastq(bam_path, out_fq, paired):
    """Reads (and mates) with any alignment on the target contig, for the control re-map."""
    names = os.path.join(W, "names.txt")
    sh(f"samtools view {bam_path} '{contig}' | cut -f1 | sort -u > {names}")
    if paired:
        sh(f"samtools view -h -N {names} {bam_path} | samtools sort -n -@4 - "
           f"| samtools fastq -F 0x900 -1 {out_fq}_1.fq -2 {out_fq}_2.fq -s /dev/null -0 /dev/null - 2>/dev/null")
    else:
        sh(f"samtools view -h -N {names} {bam_path} | samtools fastq -F 0x900 - > {out_fq}.fq 2>/dev/null")


# ---- 3./4. metrics --------------------------------------------------------------------
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def metrics(bam_path, kind, ref_seq, other_blocks=None):
    """Coverage, physical continuity, clip pile-ups and discordance over the gene.
    Depth is counted twice: all primary alignments (any MAPQ) and unique ones (MAPQ >= 20).
    A window where unique depth is < 10 % of all-MAPQ depth (and all-MAPQ depth >= 10) is a
    multi-copy (repeat) window: short reads cannot place it, so it is reported, excluded
    from the break test, and checked for long reads spanning it instead. Boundary classes and
    the break test are set afterwards by classify()."""
    lo, hi = gs, ge
    n_win = (hi - lo + 1 + WIN - 1) // WIN
    bounds = [lo + WIN * i for i in range(1, n_win)]
    dep_all = [0] * (hi - lo + 2)
    dep_q = [0] * (hi - lo + 2)
    disc = [0] * n_win
    pairs_in = [0] * n_win
    clips = {}
    inserts, frags, rlens = [], [], []
    maxdiv = MAX_DIV if kind == "illumina" else 0.15
    p = subprocess.Popen(f"samtools view -F 0x904 {bam_path} '{contig}:{max(1, lo - 20000)}-{hi + 20000}'", shell=True,
                         stdout=subprocess.PIPE, text=True, executable="/bin/bash")
    for line in p.stdout:
        f = line.split("\t")
        flag, pos, mapq, cig = int(f[1]), int(f[3]), int(f[4]), f[5]
        nm = next((int(x[5:]) for x in f[11:] if x.startswith("NM:i:")), 0)
        if len(rlens) < 20000 and not cig.count("H"):
            rlens.append(len(f[9]))
        alen, end, lclip, rclip = 0, pos, 0, 0
        ops = CIG.findall(cig)
        for i, (n, op) in enumerate(ops):
            n = int(n)
            if op in "M=X":
                alen += n
                end += n
            elif op == "D":
                end += n
            elif op == "S":
                if i == 0:
                    lclip = n
                else:
                    rclip = n
        if alen == 0 or nm / alen > maxdiv:
            continue
        uniq = mapq >= MINQ
        a, b = max(lo, pos), min(hi, end - 1)
        if a <= b:
            dep_all[a - lo] += 1
            dep_all[b - lo + 1] -= 1
            if uniq:
                dep_q[a - lo] += 1
                dep_q[b - lo + 1] -= 1
        if not uniq:
            continue
        if lclip >= 20 and lo <= pos <= hi:
            clips[pos] = clips.get(pos, 0) + 1
        if rclip >= 20 and lo <= end - 1 <= hi:
            clips[end - 1] = clips.get(end - 1, 0) + 1
        if kind == "illumina":
            if lo <= pos <= hi:
                w = (pos - lo) // WIN
                pairs_in[w] += 1
                mate_other = f[6] not in ("=", contig)
                if mate_other or not (flag & 0x2):
                    disc[w] += 1
            if (flag & 0x2) and f[6] in ("=", contig) and int(f[8]) > 0:
                inserts.append(int(f[8]))
                frags.append((pos, pos + int(f[8]) - 1))
        else:
            frags.append((pos, end - 1))
    p.wait()
    med_ins = statistics.median(inserts) if inserts else 0
    margin = 20 if kind == "illumina" else 500
    if kind == "illumina" and med_ins:
        frags = [x for x in frags if x[1] - x[0] + 1 <= 3 * med_ins]
    starts = sorted(a for a, b in frags)

    def spanning(x, m=margin):
        return sum(1 for a, b in frags if a <= x - m and b >= x + m)
    sdiff = [0] * (len(bounds) + 1)
    for a, b in frags:
        k1 = bisect_left(bounds, a + margin)
        k2 = bisect_right(bounds, b - margin)
        if k1 < k2:
            sdiff[k1] += 1
            sdiff[k2] -= 1
    span, acc = [], 0
    for k in range(len(bounds)):
        acc += sdiff[k]
        span.append(acc)

    def cum(d):
        run, out = 0, []
        for x in range(hi - lo + 1):
            run += d[x]
            out.append(run)
        return out
    da, dq = cum(dep_all), cum(dep_q)
    wa = [round(sum(da[i:i + WIN]) / len(da[i:i + WIN]), 1) for i in range(0, len(da), WIN)]
    wq = [round(sum(dq[i:i + WIN]) / len(dq[i:i + WIN]), 1) for i in range(0, len(dq), WIN)]
    rep = [wa[i] >= 10 and wq[i] < 0.1 * wa[i] for i in range(len(wa))]
    blocks, k = [], 0
    while k < len(rep):
        if rep[k]:
            j = k
            while j + 1 < len(rep) and rep[j + 1]:
                j += 1
            bs, be = lo + k * WIN, min(hi, lo + (j + 1) * WIN - 1)
            blocks.append({"start": bs, "end": be, "kb": round((be - bs + 1) / 1000, 1),
                           "all_mapq_depth": round(sum(wa[k:j + 1]) / (j - k + 1), 1)})
            k = j + 1
        else:
            k += 1
    uniq_win = [i for i in range(len(wq)) if not rep[i]]
    uq_depths = [dq[x] for i in uniq_win for x in range(i * WIN, min(len(dq), (i + 1) * WIN))]
    gene_med = statistics.median(uq_depths) if uq_depths else 0
    out = subprocess.run(f"samtools depth -a -Q {MINQ} -r '{contig}' {bam_path}", shell=True, capture_output=True,
                         text=True, executable="/bin/bash").stdout.split("\n")
    contig_depth = []
    for i, line in enumerate(out):
        if i % 50 == 0 and line:
            c, ps, d = line.split("\t")
            if not (lo - 5000 <= int(ps) <= hi + 5000):
                contig_depth.append(int(d))
    ctg_med = statistics.median(contig_depth) if contig_depth else 0
    med_span = statistics.median(span) if span else 0
    piles = []
    for p_, n in clips.items():
        d = dq[min(len(dq) - 1, max(0, p_ - lo))]
        if n >= max(5, 0.2 * d):
            sp = spanning(p_)
            piles.append({"pos": p_, "clipped_reads": n, "depth": d, "spanning_at_pos": sp,
                          "continuity_lost": sp < 0.2 * med_span})
    piles.sort(key=lambda x: -x["clipped_reads"])
    need_depth = 10 if kind == "illumina" else 5
    res = {
        "kind": kind, "gene_bp": hi - lo + 1, "windows": n_win,
        "gene_median_depth_unique_windows": gene_med, "contig_median_depth_outside_gene": ctg_med,
        "depth_ratio_gene_vs_contig": round(gene_med / ctg_med, 2) if ctg_med else None,
        "frac_gene_depth_ge5_all_mapq": round(sum(1 for d in da if d >= 5) / len(da), 4),
        "frac_unique_windows_depth_ge5": round(sum(1 for i in uniq_win for x in range(i * WIN, min(len(dq), (i + 1) * WIN))
                                                   if dq[x] >= 5) / max(1, len(uq_depths)), 4),
        "repeat_blocks": blocks, "repeat_kb": round(sum(b_["kb"] for b_ in blocks), 1),
        "window_depth_unique": wq, "window_depth_all_mapq": wa,
        "boundaries": len(bounds), "median_spanning_per_boundary": med_span,
        "spanning_per_boundary": span,
        "clip_pileups": piles[:20],
        "clip_pileups_with_continuity_lost": [x for x in piles if x["continuity_lost"]][:10],
    }
    if kind == "illumina":
        res["median_insert"] = med_ins
        res["max_discordant_frac_window"] = round(max((disc[i] / pairs_in[i]) for i in range(n_win)
                                                      if pairs_in[i] >= 20), 3) if any(x >= 20 for x in pairs_in) else None
    else:
        for b_ in blocks:      # long reads that span each repeat block with >= 500 bp unique flank each side
            b_["ont_reads_spanning"] = spanning_block(frags, b_["start"], b_["end"])
        # the short-read repeat blocks: long reads (MAPQ >= 20) spanning each with >= 500 bp flanks
        res["short_read_repeat_blocks_spanned_by_long_reads"] = [
            {"start": b_["start"], "end": b_["end"], "ont_reads_spanning": spanning_block(frags, b_["start"], b_["end"])}
            for b_ in (other_blocks or [])]
    res["_frags"] = frags
    if rlens:
        res["read_length"] = int(statistics.median(rlens))
    return res


def orf_consensus(bam_path, ref_seq):
    """CDS positions where >= 50 % of reads (depth >= 5) disagree with the reference by an
    indel (frameshift candidate) or a substitution that makes an in-frame stop."""
    cds = (gs, ge)
    out = subprocess.Popen(f"samtools mpileup -B -Q 13 -q {MINQ} -r '{contig}:{gs}-{ge}' -f {ref} {bam_path} 2>/dev/null",
                           shell=True, stdout=subprocess.PIPE, text=True, executable="/bin/bash")
    indel_major, stop_major, low = [], [], 0
    for line in out.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        pos, d, bases = int(f[1]), int(f[3]), f[4]
        if d < 5:
            low += 1
            continue
        n_indel = len(re.findall(r"[+-]\d+", bases))
        if n_indel >= 0.5 * d:
            indel_major.append(pos)
        clean = re.sub(r"\^.|\$|[+-](\d+)", "", bases)
        alt = {}
        for ch in clean.upper():
            if ch in "ACGT":
                alt[ch] = alt.get(ch, 0) + 1
        if alt:
            b, n = max(alt.items(), key=lambda x: x[1])
            if n >= 0.5 * d:
                # does this substitution create an in-frame stop on the coding strand?
                off = (pos - cds[0]) if job["strand"] == "+" else (cds[1] - pos)
                cs = off - off % 3
                if job["strand"] == "+":
                    codon = list(ref_seq[cds[0] - 1 + cs: cds[0] - 1 + cs + 3])
                    codon[off % 3] = b
                else:
                    seg = ref_seq[cds[1] - cs - 3: cds[1] - cs].translate(COMP)[::-1]
                    codon = list(seg)
                    codon[off % 3] = b.translate(COMP)
                if "".join(codon) in STOPS:
                    stop_major.append(pos)
    out.wait()
    return {"indel_majority_positions": indel_major[:20], "n_indel_majority": len(indel_major),
            "stop_creating_majority_positions": stop_major[:20], "n_stop_majority": len(stop_major),
            "cds_positions_depth_lt5": low}


result = {"job_id": JID, "contig": contig, "contig_length": clen, "gene": [gs, ge], "strand": job["strand"],
          "chimera_junction": mid, "thresholds": VERDICT, "min_mapq": MINQ, "max_divergence": MAX_DIV, "rules": RULES}
ref_seq = seqs[contig]
if args.reverdict:
    # re-judge a stored result under the current rules (no reads): the window arrays and the
    # spanning counts are stored; the library mappability is recomputed from the assembly
    old = json.load(open(args.reverdict))
    prev = old.get("previous") or {"rules": old.get("rules", "v2 (1-kb read-depth repeat mask only)"),
                                   "final_verdict": old.get("final_verdict"),
                                   "verdicts": {k: v["verdict"] for k, v in old.items()
                                                if isinstance(v, dict) and "verdict" in v}}
    result.update({k: v for k, v in old.items() if k not in result or k == "wall_s"})
    result.update({"previous": prev, "rules": RULES, "rejudged_without_reads": True})
    rl = args.read_len or job.get("illumina_read_length") or 150
    if "illumina" in result:
        judge_illumina(result["illumina"], gs, ge, ref, ref_seq, rl)
    if "negative_control_illumina" in result:
        judge_illumina(result["negative_control_illumina"], CTL[0], CTL[1], chim, CHIM_SEQ, rl)
    for k in ("ont", "negative_control_ont"):
        if k in result:
            classify(result[k], "ont", gs if k == "ont" else CTL[0])
else:
    if os.path.exists(bam + ".ok"):
        m = metrics(bam, "illumina", ref_seq)
        result["illumina"] = judge_illumina(m, gs, ge, ref, ref_seq, m.get("read_length", 150), m.pop("_frags"))
        result["illumina"]["orf"] = orf_consensus(bam, ref_seq)
        result["illumina_reads_processed"] = processed(f"{W}/mm_illumina.log")
        result["illumina_expected_reads"] = job.get("illumina_expected_reads")
        # negative control on the chimeric reference
        reads_on_contig_to_fastq(bam, f"{W}/ctl", True)
        cbam = os.path.join(W, "chimera.bam")
        sh(f"minimap2 -t {args.threads} -ax sr {chim} {W}/ctl_1.fq {W}/ctl_2.fq 2>/dev/null | samtools sort -@4 -o {cbam} - "
           f"&& samtools index {cbam}")
        gsave = (gs, ge)
        gs, ge = CTL                             # look at the junction neighbourhood
        m = metrics(cbam, "illumina", ref_seq)
        result["negative_control_illumina"] = judge_illumina(m, gs, ge, chim, CHIM_SEQ, m.get("read_length", 150),
                                                             m.pop("_frags"))
        gs, ge = gsave
    if os.path.exists(ont_bam + ".ok"):
        ill_ = result.get("illumina", {})
        unresolved = ([{"start": b["start"], "end": b["end"]} for b in ill_.get("repeat_blocks", [])]
                      + [{"start": b["start"], "end": b["end"]} for b in ill_.get("library_unplaceable_blocks", [])
                         if not b["crossed"]])
        m = metrics(ont_bam, "ont", ref_seq, other_blocks=unresolved)
        result["ont"] = classify(m, "ont", gs, None, m.pop("_frags"))
        result["ont_reads_processed"] = processed(f"{W}/mm_ont.log")
        result["ont_expected_reads"] = job.get("ont_expected_reads")
        reads_on_contig_to_fastq(ont_bam, f"{W}/ctlont", False)
        cbam = os.path.join(W, "chimera_ont.bam")
        sh(f"minimap2 -t {args.threads} -ax map-ont {chim} {W}/ctlont.fq 2>/dev/null | samtools sort -@4 -o {cbam} - "
           f"&& samtools index {cbam}")
        gsave = (gs, ge)
        gs, ge = CTL
        m = metrics(cbam, "ont", ref_seq)
        result["negative_control_ont"] = classify(m, "ont", gs, None, m.pop("_frags"))
        gs, ge = gsave


def verdict(m, kind):
    """supported / against / inconclusive (thresholds in VERDICT); a missing signal counts
    against the assembly only where coverage is adequate and the region is single-copy."""
    if not m:
        return "no data"
    need_depth = 10 if kind == "illumina" else 5
    adequate = (m["gene_median_depth_unique_windows"] or 0) >= need_depth
    why_against, why_weak, notes = [], [], []
    if m.get("break_boundaries"):
        why_against.append(f"{len(m['break_boundaries'])} single-copy boundary(ies) with 0 spanning "
                           + ("pairs" if kind == "illumina" else "reads") + " between well-covered windows")
    if adequate and m["clip_pileups_with_continuity_lost"]:
        why_against.append(f"{len(m['clip_pileups_with_continuity_lost'])} clip pile-up(s) where physical coverage drops")
    if kind == "illumina" and adequate and (m.get("max_discordant_frac_window") or 0) > VERDICT["max_discordant_frac"]:
        why_against.append(f"discordant pairs up to {m['max_discordant_frac_window']:.0%} in a window")
    r = m.get("depth_ratio_gene_vs_contig")
    if r is not None and (m.get("contig_median_depth_outside_gene") or 0) >= need_depth and \
            not (VERDICT["depth_ratio"][0] <= r <= VERDICT["depth_ratio"][1]):
        why_against.append(f"gene/contig depth ratio {r}")
    if kind == "illumina" and adequate and m.get("orf", {}).get("n_indel_majority"):
        why_against.append(f"{m['orf']['n_indel_majority']} majority-indel position(s) in the CDS")
    if kind == "illumina" and adequate and m.get("orf", {}).get("n_stop_majority"):
        why_against.append(f"{m['orf']['n_stop_majority']} stop-creating majority substitution(s)")
    if not adequate:
        why_weak.append(f"unique-window median depth {m['gene_median_depth_unique_windows']} < {need_depth}")
    if m["frac_unique_windows_depth_ge5"] < VERDICT["min_frac_depth5"]:
        why_weak.append(f"only {m['frac_unique_windows_depth_ge5']:.1%} of single-copy gene bases at depth >= 5")
    need = VERDICT["min_pairs_per_boundary"] if kind == "illumina" else VERDICT["min_ont_per_boundary"]
    mn = m.get("min_spanning_unique_boundaries")
    if mn is not None and mn < need:
        why_weak.append(f"weakest single-copy boundary spanned by {mn}")
    if m["repeat_blocks"]:
        if kind == "ont":
            unspanned = [b for b in m["repeat_blocks"] if b.get("ont_reads_spanning", 0) < VERDICT["min_ont_per_boundary"]]
            if unspanned:
                why_weak.append(f"{len(unspanned)} repeat block(s) not spanned by >= 3 long reads")
            else:
                notes.append(f"{len(m['repeat_blocks'])} repeat block(s) ({m['repeat_kb']} kb) spanned by long reads")
        else:
            notes.append(f"{m['repeat_kb']} kb in multi-copy block(s) that short reads cannot place")
    nc = [b for b in m.get("library_unplaceable_blocks", []) if not b.get("crossed")]
    if nc:
        notes.append(f"{m['unplaceable_kb_not_crossed']} kb in {m['library_unplaceable_blocks_n']} stretch(es) where "
                     "this library cannot place pairs uniquely (repeat), not crossed")
    note = (" [" + "; ".join(notes) + "]") if notes else ""
    if why_against:
        return "against: " + "; ".join(why_against) + note
    if why_weak:
        return "inconclusive: " + "; ".join(why_weak) + note
    return "supported by the reads" + note


for k in ("illumina", "ont", "negative_control_illumina", "negative_control_ont"):
    if k in result:
        result[k]["verdict"] = verdict(result[k], "illumina" if "illumina" in k else "ont")
        if k.startswith("negative_control"):
            result[k]["control_detected_the_join"] = result[k]["verdict"].startswith("against")
ill, ont = result.get("illumina"), result.get("ont")
fv = "no data"
if ill:
    fv = ill["verdict"]
    unresolved = ([(b["start"], b["end"]) for b in ill["repeat_blocks"]]
                  + [(b["start"], b["end"]) for b in ill.get("library_unplaceable_blocks", []) if not b["crossed"]])
    kb = union_kb(unresolved)
    if fv.startswith("supported") and unresolved:
        if ont:
            # a long read that spans a stored block with 500-bp flanks spans any stretch inside it
            sb = ont.get("short_read_repeat_blocks_spanned_by_long_reads", [])
            counts = [max([b["ont_reads_spanning"] for b in sb if b["start"] <= us and b["end"] >= ue], default=None)
                      for us, ue in unresolved]
            spanned = all(c is not None and c >= VERDICT["min_ont_per_boundary"] for c in counts) \
                and not ont["verdict"].startswith("against")
            fv = (f"supported by the reads (short reads; the {kb} kb repeat block(s) spanned by "
                  f"{min(counts)}+ long reads)" if spanned
                  else f"supported outside {kb} kb of repeat; long reads do not span it"
                  + (" (or it was not measured)" if any(c is None for c in counts) else ""))
        else:
            fv = f"supported outside {kb} kb of repeat (no long reads to resolve it)"
    if ont and ont["verdict"].startswith("against") and not fv.startswith("against"):
        fv = "conflicting: short reads " + ill["verdict"].split(":")[0] + ", long reads against"
result["final_verdict"] = fv
if args.reverdict:
    result["rejudge_wall_s"] = round(time.time() - t0)
else:
    result["wall_s"] = round(time.time() - t0)
json.dump(result, open(os.path.join(args.out, JID + ".json.tmp"), "w"), indent=1)
os.replace(os.path.join(args.out, JID + ".json.tmp"), os.path.join(args.out, JID + ".json"))
log("done", result["final_verdict"], {k: result[k]["verdict"] for k in result if isinstance(result[k], dict) and "verdict" in result[k]})
if args.reverdict and not args.work:
    shutil.rmtree(W, ignore_errors=True)
