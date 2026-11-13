"""Stage 2 grades: apply the RNA-seq verdicts (data/rnaseq_support.tsv, from
cloud/summarize.py) to the Stage 1 grades in data/locus_checks.tsv.

The Stage 1 grade is kept in column grade_stage1; grade is replaced by the Stage 2 grade
and grade_change_reason says why (or why not). Re-running starts again from grade_stage1.

Rules (verdict thresholds are in cloud/summarize.py and notes/RESULTS.md):
  contradicted                      -> D
  full (>= 90 % of introns)         A stays A; B -> A if Stage 1 already had S1 or S2
                                    (orthologue length or other assembly), else B;
                                    C -> B
  partial (50-90 %)                 C -> B; A and B unchanged
  weak, under-supported where covered (>= 20 % of the testable introns unsupported)
                                    A -> B, B -> C; C unchanged
  weak (coverage too low), no-expression, not tested, no runs: unchanged
  D stays D; for D5 (fusion) models the junctions around the fusion boundary are reported
  (an intron with 0 reads at a >= 10-fold step in junction reads supports the merge call)
"""
import csv
import glob
import json
import os
from statistics import median

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
path = os.path.join(D, "locus_checks.tsv")
rows = list(csv.DictReader(open(path), delimiter="\t"))
cols = list(rows[0].keys())
if "grade_stage1" not in cols:
    for r in rows:
        r["grade_stage1"] = r["grade"]
        r["grade_reason_stage1"] = r["grade_reason"]
    cols += ["grade_stage1", "grade_reason_stage1"]
rna = {r["md5"]: r for r in csv.DictReader(open(os.path.join(D, "rnaseq_support.tsv")), delimiter="\t")}
specs = {}
for p in ("specs.jsonl", "control_specs.jsonl"):
    fp = os.path.join(D, "loci", p)
    if os.path.exists(fp):
        for line in open(fp):
            s = json.loads(line)
            specs[s["md5"]] = s


def intron_totals(md5):
    tot = None
    for p in glob.glob(os.path.join(ROOT, "cloud_out", "results", "*", "*.json")):
        if p.endswith(".skip.json"):
            continue
        for L in json.load(open(p)).get("loci", []):
            if L["md5"] != md5:
                continue
            if tot is None:
                tot = {(i[0], i[1]): 0 for i in L["introns"]}
            for i in L["introns"]:
                tot[(i[0], i[1])] += i[2] + i[3]
    return tot


def fusion_boundary(r):
    """Junction reads around the gap between the first two fusion segments (aa)."""
    segs = [x.split(":")[0] for x in r["fusion_segments"].split(";") if ":" in x]
    if len(segs) < 2 or r["md5"] not in specs:
        return ""
    end1 = int(segs[0].split("-")[1])
    start2 = int(segs[1].split("-")[0])
    tot = intron_totals(r["md5"])
    if not tot:
        return ""
    s = specs[r["md5"]]
    ex = s["segments"]                       # transcription order
    introns, pos = [], 0
    for k in range(len(ex) - 1):
        pos += ex[k][1] - ex[k][0] + 1
        a, b = sorted([ex[k], ex[k + 1]])
        introns.append((pos // 3, tot.get((a[1] + 1, b[0] - 1), 0)))
    best = None
    for i, (aa, n) in enumerate(introns):
        if not (end1 - 1000 <= aa <= start2 + 1000):
            continue
        up = [x[1] for x in introns[max(0, i - 5):i]]
        down = [x[1] for x in introns[i + 1:i + 6]]
        if not up or not down:
            continue
        mu, md = median(up), median(down)
        step = max(mu, md) / max(1.0, min(mu, md))
        if best is None or (n, -step) < (best[1], -best[4]):
            best = (aa, n, mu, md, step)
    if not best:
        return ""
    aa, n, mu, md, step = best
    verdict = "supports the merge call" if n == 0 and step >= 10 else \
        "consistent with a merge but coverage too low to call" if n == 0 else "spliced across the boundary"
    return (f"fusion boundary: intron at aa {aa} has {n} junction reads (median {mu:g} in the 5 introns before, "
            f"{md:g} after; {step:.0f}-fold step): {verdict}")


changes = []
for r in rows:
    g0 = r["grade_stage1"]
    reason0 = r["grade_reason_stage1"]
    v = rna.get(r["md5"], {})
    verdict = v.get("verdict", "")
    g, why = g0, ""
    s1s2 = "S1 " in reason0 or "S2 " in reason0
    summary = ""
    if v.get("runs"):
        summary = (f"RNA-seq ({v['runs']} runs, {v['reads_unique']} unique reads): {v['introns_counted_supported']}/"
                   f"{v['n_introns']} introns supported ({v['introns_unique_support']} by unique reads, "
                   f"{v['introns_multi_only']} multi-mapped only{' - not counted, repeat model' if v['repeat_flag'] == 'yes' else ''}), "
                   f"CDS breadth {v['cds_breadth']}")
    if g0 == "D":
        why = "stays D" + (f"; {summary}" if summary else "")
        if "D5 " in reason0:
            fb = fusion_boundary(r)
            if fb:
                why += "; " + fb
    elif verdict == "contradicted" and v.get("contradicted_on_unique_junctions", "yes") != "yes":
        # pre-registered guard for the repeat-flag models (cloud/bundle_rep_rna/RULES.txt): missing
        # support at introns whose junction sequence repeats in the model can be a mapping artefact
        why = (f"stays {g0}: contradicted only at introns whose junction repeats in the model "
               f"(pre-registered guard: {v.get('contradicted_on_unique_junctions')}); {summary}")
    elif verdict == "contradicted":
        g, why = "D", f"D6 expressed but the annotated exon chain is not spliced: {summary}"
    elif verdict == "full":
        if g0 == "A":
            why = f"stays A; {summary}"
        elif g0 == "B":
            g = "A" if s1s2 else "B"
            why = (f"B -> A: transcript support added to the Stage 1 S1/S2 line; {summary}" if s1s2
                   else f"stays B (no S1/S2 line for A); {summary}")
        else:
            g, why = "B", f"C -> B: transcript support (S5); {summary}"
    elif verdict == "partial":
        if g0 == "C":
            g, why = "B", f"C -> B: transcript support for part of the model (S5 partial); {summary}"
        else:
            why = f"stays {g0}; {summary}"
    elif verdict.startswith("weak (under-supported where covered"):
        if g0 in ("A", "B"):
            g = {"A": "B", "B": "C"}[g0]
            why = (f"{g0} -> {g}: {v['testable_introns_unsupported']} of {v['testable_introns']} introns with "
                   f"well-covered flanking exons have no spliced reads; {summary}")
        else:
            why = f"stays {g0}; under-supported at adequate coverage; {summary}"
    elif verdict:
        why = f"stays {g0}; RNA-seq {verdict}" + (f"; {summary}" if summary else "")
    else:
        why = f"stays {g0}; not in the RNA-seq stage" + (" (no expression data for bacteria/archaea/MAGs)"
                                                        if r["superkingdom"] != "Eukaryota" else "")
    r["grade"] = g
    r["grade_stage2_rna"] = g            # the RNA-stage grade; scripts/20_apply_dna.py starts from it
    r["rnaseq_verdict"] = verdict or ("not applicable" if r["superkingdom"] != "Eukaryota" else "not selected")
    r["rnaseq_summary"] = summary
    r["grade_change_reason"] = why
    r["grade_change_reason_stage2_rna"] = why   # scripts/20_apply_dna.py starts from it
    if g != g0:
        changes.append((r["rep_accession"], g0, g))
    r["grade_reason"] = reason0 if g == g0 else f"{reason0}; Stage 2: {why}"

for c in ("rnaseq_verdict", "rnaseq_summary", "grade_change_reason", "grade_stage2_rna",
          "grade_change_reason_stage2_rna"):
    if c not in cols:
        cols.append(c)
with open(path, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")


def table(sel, title):
    grades = "ABCD"
    print(title)
    print("stage1\\stage2 " + " ".join(grades) + "  total")
    for g0 in grades:
        n = [sum(1 for r in sel if r["grade_stage1"] == g0 and r["grade"] == g) for g in grades]
        print(f"   {g0}         " + " ".join(str(x) for x in n) + f"  {sum(n)}")
    print("   total     " + " ".join(str(sum(1 for r in sel if r["grade"] == g)) for g in grades))


main = [r for r in rows if not r["md5"].startswith("control_")]
table(main, "All 510 sequences")
table([r for r in main if int(r["length"]) > 45212], "Sequences > 45,212 aa")
print("changes:", len(changes))
for c in changes:
    print("  ", c)
