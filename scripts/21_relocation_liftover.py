"""Lift-over for RNA-stage loci that were relocated between non-identical sequences.

cloud/prep_genome.py places a locus on the downloaded assembly by one 300-bp probe and moves
every exon by the same offset. That is exact when the two sequences are the same (a WGS
scaffold and its CM chromosome, a GenBank chromosome and its RefSeq copy), but not across
haplotypes: CAM4536479.1 (Lepidochelys kempii, hap2 contig OZ223597.1) was placed on the hap1
chromosome NC_136390.1, its intron coordinates drifted by 1-5 bp past the first indel between
the haplotypes, and exact-match junction counting found 13 of 334 introns (verdict
"contradicted"), although the reads carry the same junctions a few bp away.

For every relocation in the RNA driver logs, this script aligns the model's source locus to
the target locus (local copies in data/loci/dna, minimap2 -x asm5 --cs), lifts each intron
end, and compares with the coordinates the box used (constant offset). Models whose
coordinates change are written to data/followups/relocation_liftover.json; cloud/summarize.py
then counts junction reads at the lifted coordinates (annotated junctions of every model of
that species and run, plus the run's novel junctions; unique reads).
"""
import glob
import json
import os
import re
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
LOGS = sorted(glob.glob(os.path.join(ROOT, "cloud_out", "driver*.log")))
DNA = os.path.join(ROOT, "data", "loci", "dna")
OUT = os.path.join(ROOT, "data", "followups", "relocation_liftover.json")

reloc = {}                                   # species key -> {src: dst}
for p in LOGS:
    for line in open(p):
        m = re.search(r"prep ready (\S+) \d+s (\{.*\})", line)
        if m:
            for x in json.loads(m.group(2)).get("relocated", []):
                src, dst = x.split("->")
                reloc.setdefault(m.group(1), {})[src] = dst


def locus_file(contig, lo, hi):
    for p in glob.glob(os.path.join(DNA, contig + "_*.fa.gz")):
        a, b = map(int, os.path.basename(p)[len(contig) + 1:-6].split("_"))
        if a <= lo and b >= hi:
            return p, a
    return None, None


def lift_map(src_fa, dst_fa):
    """source position (1-based, in the source file's contig coordinates) -> target position."""
    out = subprocess.run(["minimap2", "-x", "asm5", "-c", "--cs", dst_fa[0], src_fa[0]],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    out = [l.split("\t") for l in out]
    if not out:
        return None, None
    f = max(out, key=lambda x: int(x[9]))
    if f[4] != "+":
        return None, "reverse-strand alignment (not handled)"
    qs, ts = int(f[2]), int(f[7])
    cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
    m, q, t = {}, qs, ts
    for op in re.findall(r"(:\d+|\*[a-z]{2}|\+[a-z]+|-[a-z]+)", cs):
        if op[0] == ":":
            n = int(op[1:])
            for i in range(n):
                m[src_fa[1] + q + i] = dst_fa[1] + t + i
            q, t = q + n, t + n
        elif op[0] == "*":
            m[src_fa[1] + q] = dst_fa[1] + t
            q, t = q + 1, t + 1
        elif op[0] == "+":
            q += len(op) - 1
        else:
            t += len(op) - 1
    return m, f"{int(f[9])}/{int(f[10])} identical bases"


result, report = {}, []
for sp, moves in sorted(reloc.items()):
    loci = []
    for b in ("bundle", "bundle_followup_rna"):
        p = os.path.join(ROOT, "cloud", b, "loci", sp + ".json")
        if os.path.exists(p):
            loci += json.load(open(p))
    res = sorted(glob.glob(os.path.join(ROOT, "cloud_out", "results", sp, "*.json")))
    if not res:
        continue
    for L in loci:
        if L["contig"] not in moves:
            continue
        dst = moves[L["contig"]]
        lo, hi = min(min(s) for s in L["segments"]), max(max(s) for s in L["segments"])
        box = next((x for f in res for x in json.load(open(f)).get("loci", []) if x["md5"] == L["md5"]), None)
        if not box:
            continue
        box_lo = min(e[0] for e in box["exons"])
        src_p, src_a = locus_file(L["contig"], lo, hi)
        span = (box_lo, box_lo + (hi - lo))
        dst_p, dst_a = locus_file(dst, *span)
        row = {"species_key": sp, "rep_accession": L["rep_accession"], "relocation": f"{L['contig']}->{dst}"}
        if not src_p or not dst_p:
            row["status"] = "no local sequence for the source or target locus: not checked"
            report.append(row)
            continue
        m, how = lift_map((src_p, src_a), (dst_p, dst_a))
        if m is None:
            row["status"] = f"no alignment ({how})"
            report.append(row)
            continue
        # the box applied one offset; a strand flip would show as a reverse alignment (refused above)
        off = box_lo - lo
        lifted, changed, missing = [], 0, 0
        for a, b in [(i[0], i[1]) for i in box["introns"]]:
            la, lb = m.get(a - off), m.get(b - off)
            if la is None or lb is None:
                missing += 1
                lifted.append(None)
                continue
            lifted.append([la, lb])
            changed += (la, lb) != (a, b)
        row.update({"alignment": how, "introns": len(lifted), "introns_moved_by_liftover": changed,
                    "introns_not_liftable": missing})
        row["status"] = "coordinates change: lift-over applied" if changed else "identical placement: no change"
        report.append(row)
        if changed:
            result[L["md5"]] = {"rep_accession": L["rep_accession"], "species_key": sp,
                                "relocation": row["relocation"], "alignment": how,
                                "introns_lifted": lifted, "introns_moved": changed}

os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(result, open(OUT, "w"), indent=0)
for r in report:
    print(r)
print(f"{len(result)} model(s) with lifted coordinates -> {OUT}")
