"""Follow-up (b) and (c): RNA-seq job bundle for
  (b) the three models lost to the job-key collision, and
  (c) the four titin models whose runs stopped on the saturation cap (Rana, Microcaecilia,
      Pantherophis, Solea): only the runs that stopped on saturation, re-run with the
      saturation rule switched off (followup_box.sh passes --sat-cov 0), so the whole run
      (up to the 60M-reads-per-file and 30-min caps) is read.
(b) details: the three models were lost to the job-key collision
(two job groups of one species got the same species key, so the second overwrote the first
group's loci file). Same genome and the same runs as the original species job; new species
key "<key>_extra" so the RNA driver does not skip it as already done.
Output: cloud/bundle_followup_rna/{jobs.tsv, loci/, exons/}
"""
import csv
import gzip
import json
import os
import shutil

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "cloud", "bundle_followup_rna")
LOST = {"XP_084498969.1": "Lepidochelys_olivacea_965140295.1",
        "CAM4536708.1": "Lepidochelys_kempii_965140285.1",
        "CAM4536479.1": "Lepidochelys_kempii_965140285.1"}
shutil.rmtree(OUT, ignore_errors=True)
os.makedirs(os.path.join(OUT, "loci"))
os.makedirs(os.path.join(OUT, "exons"))
jobs = {r["species_key"]: r for r in csv.DictReader(open(os.path.join(ROOT, "cloud", "bundle", "jobs.tsv")),
                                                     delimiter="\t")}
specs = {}
for line in open(os.path.join(D, "loci", "specs.jsonl")):
    s = json.loads(line)
    specs[s["rep_accession"]] = s
idx = json.load(open(os.path.join(D, "loci", "dna_index.json")))
checks = {r["rep_accession"]: r for r in csv.DictReader(open(os.path.join(D, "locus_checks.tsv")), delimiter="\t")}


def region(s):
    for lid, v in idx.items():
        if s["md5"] in v["md5s"]:
            seq = "".join(l.strip() for l in gzip.open(os.path.join(D, "loci", "dna", lid + ".fa.gz"), "rt")
                          if not l.startswith(">"))
            return seq, v["start"]
    raise SystemExit("no DNA for " + s["rep_accession"])


by_sp = {}
for acc, sp in LOST.items():
    by_sp.setdefault(sp, []).append(acc)
rows = []
for sp, accs in by_sp.items():
    key = sp + "_extra"
    loci = []
    with open(os.path.join(OUT, "exons", key + ".fa"), "w") as fe:
        for acc in accs:
            s = specs[acc]
            seq, start = region(s)
            lo = min(a for a, b in s["segments"])
            loci.append({"locus_id": s["md5"][:12], "md5": s["md5"], "rep_accession": acc, "contig": s["contig"],
                         "strand": s["strand"], "segments": s["segments"], "grade0": checks[acc]["grade"],
                         "probe": seq[lo - start:lo - start + 300], "probe_pos": lo})
            iv = sorted((max(1, a - 100), b + 100) for a, b in s["segments"])
            merged = []
            for a, b in iv:
                if merged and a <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            for a, b in merged:
                fe.write(f">{s['md5'][:12]}_{a}\n{seq[max(0, a - start):b - start + 1]}\n")
    json.dump(loci, open(os.path.join(OUT, "loci", key + ".json"), "w"))
    j = dict(jobs[sp])
    j["species_key"] = key
    j["role"] = "collision-rerun"
    j["n_loci"] = str(len(loci))
    j["max_locus_aa"] = str(max(int(checks[a]["length"]) for a in accs))
    rows.append(j)
    print(key, accs, j["genome_url"].split("/")[-1], json.loads(j["runs_json"]).__len__(), "runs",
          j["runs_gb"], "GB")
# (c) saturation-limited titin runs
NOSAT = {"XP_040214433.1": ["SRR11300562", "SRR7140896", "SRR7251162"],
         "XP_030066708.1": ["SRR5591423"],
         "XP_087738184.1": ["SRR32343001", "SRR9596706"],
         "KAG7502825.1": ["SRR18693561", "SRR18693670"]}
orig_loci = {}
for sp in jobs:
    p = os.path.join(ROOT, "cloud", "bundle", "loci", sp + ".json")
    if os.path.exists(p):
        for L in json.load(open(p)):
            orig_loci[L["rep_accession"]] = (sp, L)
for acc, runs in NOSAT.items():
    sp, L = orig_loci[acc]
    key = sp + "_nosat"
    json.dump([L], open(os.path.join(OUT, "loci", key + ".json"), "w"))
    shutil.copy(os.path.join(ROOT, "cloud", "bundle", "exons", sp + ".fa"), os.path.join(OUT, "exons", key + ".fa"))
    j = dict(jobs[sp])
    sel = [r for r in json.loads(j["runs_json"]) if r["run"] in runs]
    j.update({"species_key": key, "role": "no-saturation rerun", "n_loci": "1",
              "max_locus_aa": checks[acc]["length"], "runs_json": json.dumps(sel, separators=(",", ":")),
              "n_runs": str(len(sel)), "runs_gb": str(round(sum(r["bytes"] for r in sel) / 1e9, 2))})
    rows.append(j)
    print(key, acc, len(sel), "runs", j["runs_gb"], "GB")
cols = list(rows[0].keys())
with open(os.path.join(OUT, "jobs.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
