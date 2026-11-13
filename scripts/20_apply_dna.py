"""Apply the DNA re-mapping checks (data/dna_checks/*.json, cloud/dna_check.py) to the grades
of the bacterial and archaeal giants in data/locus_checks.tsv.

Section 5 of notes/FEASIBILITY.md defines grade A for bacteria and MAGs through the reads:
no N, coding-frame signal, reads spanning every window with no break or frame-shifting
indel, plus one independent line (orthologue of consistent length or another assembly).
Rules here (applied after scripts/15_apply_rnaseq.py; the Stage 2 RNA grade is kept in
grade_stage2_rna):
  reads "supported" (repeat blocks, if any, spanned by long reads) + S1 or S2 -> A
  reads "supported" without S1/S2                                           -> B (reads are the support line)
  reads "against"                                                             -> D (D7 assembly not supported by its reads)
  "supported outside N kb of repeat" or "inconclusive"                        -> unchanged, noted
Each result must be judged under the current rules of cloud/dna_check.py (rules v3 since
12 Nov: the break test only where the library can place pairs; results from older rules are
re-judged with dna_check.py --reverdict before they are copied to data/dna_checks/). A
result whose "rules" field is missing or older is refused.
The job -> sequence mapping is by the CDS coordinates on the contig.
"""
import csv
import glob
import json
import os

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
path = os.path.join(D, "locus_checks.tsv")
rows = list(csv.DictReader(open(path), delimiter="\t"))
cols = list(rows[0].keys())
# grade_stage2_rna is written by scripts/15_apply_rnaseq.py (the RNA-stage grade); this script
# always starts from it, so re-running it, or re-running 15 first, never stacks the DNA rules
if "grade_stage2_rna" not in cols or "grade_change_reason_stage2_rna" not in cols:
    raise SystemExit("run scripts/15_apply_rnaseq.py first (it writes grade_stage2_rna and its reason)")
for c in ("dna_check_verdict", "dna_check_job"):
    if c not in cols:
        cols.append(c)
specs = {}
for line in open(os.path.join(D, "loci", "specs.jsonl")):
    s = json.loads(line)
    if s["segments"]:
        specs[s["md5"]] = s
results = {}
for p in glob.glob(os.path.join(D, "dna_checks", "*.json")):
    d = json.load(open(p))
    if not str(d.get("rules", "")).startswith("v3"):
        raise SystemExit(f"{p}: judged under old rules ({d.get('rules', 'v2')}); re-judge with dna_check.py --reverdict")
    results[(d["contig"].split(".")[0], d["gene"][0], d["gene"][1])] = d
# jobs sent or run without a result yet (shown as pending)
jobs = {}
for p in glob.glob(os.path.join(ROOT, "cloud", "bundle_dna", "jobs", "*.json")) + \
        glob.glob(os.path.join(ROOT, "cloud", "bundle_dna", "done_on_mac", "*.json")):
    j = json.load(open(p))
    jobs[(j["contig"].split(".")[0], min(j["gene_start"], j["gene_end"]), max(j["gene_start"], j["gene_end"]))] = j["job_id"]
changes = []
for r in rows:
    g = r["grade_stage2_rna"]
    s = specs.get(r["md5"])
    r["grade"] = g
    r["dna_check_verdict"] = ""
    r["dna_check_job"] = ""
    r["grade_change_reason"] = r.get("grade_change_reason_stage2_rna", r.get("grade_change_reason", ""))
    if not s or len(s["segments"]) != 1:
        continue
    a, b = s["segments"][0]
    key = (s["contig"].split(".")[0], min(a, b), max(a, b))
    d = results.get(key)
    if not d:
        if key in jobs:
            r["dna_check_verdict"] = "pending: job sent, no result yet"
            r["dna_check_job"] = jobs[key]
        continue
    fv = d["final_verdict"]
    r["dna_check_verdict"] = fv
    r["dna_check_job"] = d["job_id"]
    reason = r["grade_reason_stage1"]
    s12 = "S1 " in reason or "S2 " in reason
    ill = d.get("illumina", {})
    detail = (f"reads: {d.get('illumina_reads_processed')} Illumina reads of {d.get('illumina_expected_reads')}, "
              f"median depth {ill.get('gene_median_depth_unique_windows')} (contig {ill.get('contig_median_depth_outside_gene')}), "
              f"min spanning pairs {ill.get('min_spanning_unique_boundaries')}, no break, no frame-shifting indel"
              + (f"; long reads {d.get('ont_reads_processed')}" if d.get("ont") else ""))
    if fv.startswith("supported by the reads"):
        if g in ("B", "C"):
            new = "A" if s12 else "B"
            r["grade"] = new
            r["grade_change_reason"] = (f"{g} -> {new}: assembly supported by its reads (section 5 bacterial criterion 3)"
                                        f"{' plus the Stage 1 S1/S2 line' if s12 else ''}; {fv}; {detail}")
    elif fv.startswith("against"):
        r["grade"] = "D"
        r["grade_change_reason"] = f"{g} -> D: D7 assembly not supported by its reads: {fv}"
    else:
        r["grade_change_reason"] = (r.get("grade_change_reason", "") + f"; reads: {fv}").strip("; ")
    if r["grade"] != g:
        changes.append((r["rep_accession"], r["organism"], r["length"], g, r["grade"]))
    r["grade_reason"] = r["grade_reason_stage1"] + "; " + r["grade_change_reason"]
with open(path, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in cols) + "\n")
print("DNA checks applied:", len(results), "jobs; grade changes:", changes)
