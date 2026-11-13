"""Threshold sensitivity of the grade counts (Stage 1 and RNA-seq rules).

The real pipeline (scripts/13_grade.py -> cloud/summarize.py -> scripts/15_apply_rnaseq.py ->
scripts/20_apply_dna.py) is run in a throw-away copy of the repository tree, once with the
published thresholds (must reproduce data/locus_checks.tsv) and once per varied threshold,
one at a time. Stage 1 thresholds are varied through environment overrides of
scripts/13_grade.py, RNA thresholds through cloud/summarize.py's command-line options.
The DNA-check thresholds are varied separately (controls/dna/sensitivity.tsv).
Output: controls/sensitivity/sensitivity_stage12.tsv
"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
PKZ = 45212

VARIANTS = [("published thresholds", {}, []),
            ("D1 rule v1: repeat pattern alone is a contradiction (used until 13 Nov)", {"GP_D1_RULE": "v1"}, []),
            ("RNA verdict order v1: covered-intron test only below 50 % support (used until 13 Nov)", {},
             ["--verdict-order", "v1"])]
for v in ("0.5", "0.7", "0.8"):
    VARIANTS.append((f"D1 period identity >= {v} (published 0.6)", {"GP_D1_ID": v}, []))
    VARIANTS.append((f"D1 period identity >= {v}, D1 rule v1", {"GP_D1_ID": v, "GP_D1_RULE": "v1"}, []))
for v in ("0.6", "0.9"):
    VARIANTS.append((f"D1 rep25 >= {v} (published 0.8)", {"GP_D1_REP25": v}, []))
    VARIANTS.append((f"D1 rep25 >= {v}, D1 rule v1", {"GP_D1_REP25": v, "GP_D1_RULE": "v1"}, []))
for v in ("0.02", "0.10", "0.20"):
    VARIANTS.append((f"D1-D3 domain guard < {v} (published 0.05)", {"GP_DOM": v}, []))
for v in ("0.3", "0.7"):
    VARIANTS.append((f"D2 non-coding windows >= {v} (published 0.5)", {"GP_D2_NC": v}, []))
for v in ("1", "3"):
    VARIANTS.append((f"D3 open frames >= {v} (published 2)", {"GP_D3_FRAMES": v}, []))
for v in ("0.90", "1.00"):
    VARIANTS.append((f"coding signal Fickett median >= {v} (published 0.95)", {"GP_FICKETT_MED": v}, []))
VARIANTS.append(("record twins counted as independent (S1 species key and S2 assembly key before 13 Nov)", {"GP_TWIN_FIX": "0"}, []))
VARIANTS.append(("S6 peptides: >= 3 of 5 length fifths (published: all 5; first draft 13 Nov)", {"GP_S6_FIFTHS": "3"}, []))
for v in ("1", "2", "5", "10"):
    VARIANTS.append((f"spliced reads per intron >= {v} (published 3)", {}, ["--min-reads", v]))
for v in ("0.8", "0.95"):
    VARIANTS.append((f"'full' intron fraction >= {v} (published 0.9)", {}, ["--full-frac", v]))
for v in ("0.3", "0.7"):
    VARIANTS.append((f"'partial' intron fraction >= {v} (published 0.5)", {}, ["--partial-frac", v]))
for v in ("5", "20"):
    VARIANTS.append((f"testable-intron flank depth >= {v} (published 10)", {}, ["--flank-depth", v]))


def shadow():
    t = tempfile.mkdtemp(prefix="gp_sens_")
    shutil.copytree(os.path.join(ROOT, "scripts"), os.path.join(t, "scripts"), ignore=shutil.ignore_patterns("__pycache__"))
    os.makedirs(os.path.join(t, "cloud"))
    shutil.copy(os.path.join(ROOT, "cloud", "summarize.py"), os.path.join(t, "cloud", "summarize.py"))
    for d in ("bundle", "bundle_followup_rna", "bundle_dna", "bundle_rep_rna"):
        os.symlink(os.path.join(ROOT, "cloud", d), os.path.join(t, "cloud", d))
    os.makedirs(os.path.join(t, "controls"))
    os.symlink(os.path.join(ROOT, "controls", "rna_rep"), os.path.join(t, "controls", "rna_rep"))
    os.symlink(os.path.join(ROOT, "cloud_out"), os.path.join(t, "cloud_out"))
    os.makedirs(os.path.join(t, "data"))
    for name in os.listdir(os.path.join(ROOT, "data")):
        src = os.path.join(ROOT, "data", name)
        dst = os.path.join(t, "data", name)
        if name in ("loci", "fasta", "raw", "refdb", "manual", "dna_checks", "followups", "hmm"):
            os.symlink(src, dst)
        elif os.path.isfile(src):
            shutil.copy(src, dst)
    return t


def run(t, env, rna_args):
    e = dict(os.environ, **env)
    q = dict(stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True, env=e)
    subprocess.run([PY, os.path.join(t, "scripts", "13_grade.py")], **q)
    subprocess.run([PY, os.path.join(t, "cloud", "summarize.py"), "--checks", os.path.join(t, "data", "locus_checks.tsv"),
                    "--out", os.path.join(t, "data", "rnaseq_support.tsv")] + rna_args, **q)
    subprocess.run([PY, os.path.join(t, "scripts", "15_apply_rnaseq.py")], **q)
    subprocess.run([PY, os.path.join(t, "scripts", "20_apply_dna.py")], **q)
    rows = [r for r in csv.DictReader(open(os.path.join(t, "data", "locus_checks.tsv")), delimiter="\t")
            if not r["md5"].startswith("control_")]
    return rows


def counts(rows, col):
    a = Counter(r[col] for r in rows)
    b = Counter(r[col] for r in rows if int(r["length"]) > PKZ)
    return a, b


base_real = {r["md5"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "locus_checks.tsv")), delimiter="\t")}
out = []
t = shadow()
try:
    for label, env, rna in VARIANTS:
        rows = run(t, env, rna)
        if not env and not rna:
            diff = sum(1 for r in rows if r["grade"] != base_real[r["md5"]]["grade"] or
                       r["grade_stage1"] != base_real[r["md5"]]["grade_stage1"])
            print("baseline reproduces data/locus_checks.tsv grades:", diff == 0, "(differences:", diff, ")", flush=True)
            assert diff == 0
        fa, fb = counts(rows, "grade")
        sa, sb = counts(rows, "grade_stage1")
        moved = sum(1 for r in rows if r["grade"] != base_real[r["md5"]]["grade"])
        moved_big = sum(1 for r in rows if r["grade"] != base_real[r["md5"]]["grade"] and int(r["length"]) > PKZ)
        cand = sorted(r["rep_accession"] for r in rows if int(r["length"]) > PKZ and r["grade"] in "AB")
        out.append({"variant": label,
                    **{f"all_{g}": fa[g] for g in "ABCD"}, **{f"big_{g}": fb[g] for g in "ABCD"},
                    "sequences_changing_grade": moved, "of_which_above_45212": moved_big,
                    "big_A_or_B": ";".join(cand)})
        print(label, dict(fa), dict(fb), moved, moved_big, flush=True)
finally:
    shutil.rmtree(t, ignore_errors=True)
cols = list(out[0].keys())
with open(os.path.join(OUT, "sensitivity_stage12.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in out:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
print("written", os.path.relpath(os.path.join(OUT, "sensitivity_stage12.tsv"), ROOT))
