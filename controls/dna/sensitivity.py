"""DNA-check threshold sensitivity: re-judge the five real source-read checks (data/dna_checks/*.json)
with cloud/dna_check.py --reverdict, varying one VERDICT threshold at a time from the defaults.
The implied grade follows scripts/20_apply_dna.py: supported -> A with S1/S2 (else B), against -> D,
anything else -> grade unchanged (the RNA-stage grade). Output: controls/dna/sensitivity.tsv.
Run: python3 controls/dna/sensitivity.py   (needs the assemblies in data/loci/{mags_dl,nocardia_dl})
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
D = os.path.join(ROOT, "data")
JOBS = {
    "omnitrophota_DAOOUB_HOW35280": ("cloud/bundle_dna/jobs", "data/loci/mags_dl/DAOOUB01.fasta.gz", 151),
    "omnitrophota_DAOPMG_HPD84636": ("cloud/bundle_dna/jobs", "data/loci/mags_dl/DAOPMG01.fasta.gz", 151),
    "omnitrophota_DAOTDT_HRZ03494": ("cloud/bundle_dna/jobs", "data/loci/mags_dl/DAOTDT01.fasta.gz", 151),
    "nocardia_NBC_01730_WSG58785": ("cloud/bundle_dna/done_on_mac", "data/loci/nocardia_dl/GCA_035920445.1.fna.gz", 150),
    "nocardia_NBC_00881_WSY60125": ("cloud/bundle_dna/done_on_mac", "data/loci/nocardia_dl/GCA_045207455.1.fna.gz", 150),
}
SETTINGS = [("default", {})]
for k, vals in (("min_pairs_per_boundary", (3, 10)), ("min_unique_depth_illumina", (5, 20)),
                ("min_frac_depth5", (0.95,)), ("max_discordant_frac", (0.05, 0.20))):
    for v in vals:
        SETTINGS.append((f"{k}={v}", {k: v}))

# sequence -> RNA-stage grade and S1/S2 (read-only use of data/locus_checks.tsv)
lc = {r["dna_check_job"]: r for r in csv.DictReader(open(os.path.join(D, "locus_checks.tsv")), delimiter="\t")
      if r.get("dna_check_job")}


def grade(fv, row):
    g = row["grade_stage2_rna"]
    s12 = "S1 " in row["grade_reason_stage1"] or "S2 " in row["grade_reason_stage1"]
    if fv.startswith("supported by the reads"):
        return ("A" if s12 else "B") if g in ("B", "C") else g
    if fv.startswith("against"):
        return "D"
    return g


out = tempfile.mkdtemp(prefix="dnasens_")
rows = []
for name, ov in SETTINGS:
    env = dict(os.environ)
    if ov:
        env["DNA_VERDICT_OVERRIDE"] = json.dumps(ov)
    for jid, (jdir, asm, rl) in JOBS.items():
        r = subprocess.run([sys.executable, os.path.join(ROOT, "cloud", "dna_check.py"), os.path.join(ROOT, jdir, jid + ".json"),
                            out, "--reverdict", os.path.join(D, "dna_checks", jid + ".json"),
                            "--assembly", os.path.join(ROOT, asm), "--read-len", str(rl), "--threads", "4"],
                           env=env, capture_output=True, text=True)
        if r.returncode:
            print(r.stderr[-2000:])
            raise SystemExit(f"{jid} {name} failed")
        d = json.load(open(os.path.join(out, jid + ".json")))
        row = lc[jid]
        rows.append({"setting": name, "job": jid, "accession": row["rep_accession"],
                     "final_verdict": d["final_verdict"], "grade": grade(d["final_verdict"], row),
                     "grade_now": row["grade"], "rna_stage_grade": row["grade_stage2_rna"]})
        print(name, jid, "|", d["final_verdict"][:90], "|", rows[-1]["grade"], flush=True)
with open(os.path.join(ROOT, "controls", "dna", "sensitivity.tsv"), "w") as f:
    cols = list(rows[0].keys())
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
print("written controls/dna/sensitivity.tsv")
