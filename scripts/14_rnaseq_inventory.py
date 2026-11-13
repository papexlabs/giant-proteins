"""Completeness and honesty check of the AWS RNA-seq stage (results synced to cloud_out/).

For every run in the job lists (cloud/bundle/jobs.tsv, and cloud/bundle_followup_rna/jobs.tsv
for the 12 Nov follow-up: _nosat = re-run without the saturation stop, _extra = loci lost to
the job-key collision; cloud/bundle_rep_rna for the 13 Nov repeat-flag stage, _rep; formerly:
the job-key collision) plus any extra run that has a
result (runs from the first, v1 box): status, driver version, stop rule, reads read vs
reads in the run, candidate reads, and leftovers of interrupted writes (.part/.tmp files
left in cloud_out because rsync does not delete). A run counts as "usable" if it has a
final JSON with the per-locus summary. Duplicate results (one run in two places) and
results without a job entry are reported.
Output: data/rnaseq_runs.tsv and a printed summary.
"""
import csv
import glob
import json
import os
from collections import Counter, defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..")
RES = os.path.join(ROOT, "cloud_out", "results")
jobs = []
for b in ("bundle", "bundle_followup_rna", "bundle_rep_rna"):
    jp = os.path.join(ROOT, "cloud", b, "jobs.tsv")
    if os.path.exists(jp):
        jobs += list(csv.DictReader(open(jp), delimiter="\t"))

planned = {}
for j in jobs:
    for r in json.loads(j["runs_json"]):
        planned[(j["species_key"], r["run"])] = r

files = [f for f in glob.glob(os.path.join(RES, "*", "*")) if not f.split(os.sep)[-2].startswith("_")]
results = defaultdict(list)
leftovers = []
for f in files:
    sp, name = f.split(os.sep)[-2], os.path.basename(f)
    if name.endswith((".part", ".tmp")):
        leftovers.append((sp, name))
    elif name.endswith(".skip.json"):
        results[(sp, name[:-10])].append(("skip", f))
    elif name.endswith(".json"):
        results[(sp, name[:-5])].append(("json", f))

rows = []
for key in sorted(set(planned) | set(results)):
    sp, run = key
    items = results.get(key, [])
    meta = planned.get(key, {})
    row = {"species_key": sp, "run": run, "in_job_list": "yes" if key in planned else "no (v1 box)",
           "status": "missing", "driver": "", "stop_rule": "", "reads_processed": "", "run_read_count": "",
           "frac_of_run_read": "", "candidate_reads": "", "stage_a_s": "", "skip_reason": "",
           "n_result_files": len(items)}
    for kind, f in items:
        if kind == "skip":
            row["status"] = "skipped"
            row["skip_reason"] = json.load(open(f)).get("reason", "")
            continue
        d = json.load(open(f))
        row["status"] = "usable"
        row["driver"] = "v2" if "stop_rule" in d else "v1"
        rules = d.get("stop_rule") or []
        row["stop_rule"] = "/".join(sorted(set(rules)))
        rm = d.get("run_meta", {}) or meta
        row["run_read_count"] = rm.get("read_count", "")
        if "reads_processed" in d:
            row["reads_processed"] = d["reads_processed"]
            n_files = max(1, len(meta.get("fastq_urls", []) or [1]))
            # ENA read_count counts spots (pairs); reads_processed counts reads over all files
            denom = (rm.get("read_count") or 0) * n_files
            row["frac_of_run_read"] = round(d["reads_processed"] / denom, 3) if denom else ""
        row["candidate_reads"] = d.get("candidate_reads", "")
        row["stage_a_s"] = d.get("stage_a_s", "")
    rows.append(row)

cols = list(rows[0].keys())
with open(os.path.join(ROOT, "data", "rnaseq_runs.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")

usable = [r for r in rows if r["status"] == "usable"]
print(f"planned runs {len(planned)} in {len(jobs)} species; result keys {len(results)}")
print("status:", Counter(r["status"] for r in rows))
print("usable by driver:", Counter(r["driver"] for r in usable))
print("usable not in the current job list (v1 box):", [(r["species_key"][:20], r["run"]) for r in usable
                                                       if r["in_job_list"] != "yes"])
print("missing (planned, no result):", [(r["species_key"][:25], r["run"]) for r in rows if r["status"] == "missing"])
print("skipped:", [(r["species_key"][:25], r["run"], r["skip_reason"]) for r in rows if r["status"] == "skipped"])
print("duplicate result files for one run:", [(r["species_key"], r["run"]) for r in rows if r["n_result_files"] > 1])
rule_parts = Counter()
for r in usable:
    for p in (r["stop_rule"] or "v1-no-rule").split("/"):
        rule_parts[p] += 1
print("stop-rule components among usable runs (a run can have one per FASTQ file):", rule_parts)
print("usable runs by combined stop rule:", Counter(r["stop_rule"] or "v1 (no rule recorded)" for r in usable))
fr = [r["frac_of_run_read"] for r in usable if r["frac_of_run_read"] != ""]
print(f"fraction of each run actually read (v2): median {sorted(fr)[len(fr)//2]:.2f}; "
      f"< 0.5: {sum(1 for x in fr if x < 0.5)}; >= 0.95: {sum(1 for x in fr if x >= 0.95)} of {len(fr)}")
sp_usable = Counter(r["species_key"] for r in usable)
print("species in job list without a usable run:", [j["species_key"] for j in jobs if sp_usable[j["species_key"]] == 0])
done = {os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(RES, "*", "DONE"))}
print("species without DONE:", [j["species_key"] for j in jobs if j["species_key"] not in done])
print("leftover partial files:", leftovers)
for sp, name in leftovers:
    run = name.split(".")[0]
    print(f"  {sp}/{name}: final result for the run exists: {any(k == (sp, run) for k in results)}")
