#!/usr/bin/env python3
"""Resumable driver for the DNA re-mapping checks (cloud/dna_check.py), box or Mac.

For every job in <bundle_dna>/jobs/*.json without a result in <results>/_dna/<job_id>.json,
run dna_check.py; up to --par jobs at once. Each finished job is one durable JSON (atomic
rename), which the Mac sync picks up; a relaunched box skips finished jobs. Appends to
<work>/driver.log like the RNA driver.
Usage: dna_driver.py <bundle_dna> <work> [--par 3] [--threads 32] [--max-pairs N]
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("bundle")
ap.add_argument("work")
ap.add_argument("--par", type=int, default=3)
ap.add_argument("--threads", type=int, default=os.cpu_count())
ap.add_argument("--max-pairs", type=int, default=0)
ap.add_argument("--only", default="", help="comma-separated job ids")
args = ap.parse_args()
RES = os.path.join(args.work, "results", "_dna")
TMP = os.path.join(args.work, "tmp_dna")
os.makedirs(RES, exist_ok=True)
os.makedirs(TMP, exist_ok=True)
LOG = open(os.path.join(args.work, "driver.log"), "a")


def log(*a):
    msg = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " [dna] " + " ".join(str(x) for x in a)
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()


jobs = sorted(glob.glob(os.path.join(args.bundle, "jobs", "*.json")))
if args.only:
    keep = set(args.only.split(","))
    jobs = [j for j in jobs if os.path.basename(j)[:-5] in keep]
todo = [j for j in jobs if not os.path.exists(os.path.join(RES, os.path.basename(j)))]
log(f"dna driver start: {len(jobs)} jobs, {len(todo)} to do, par {args.par}")
tpj = max(2, args.threads // max(1, args.par))
running = []
last_beat = time.time()
while todo or running:
    if time.time() - last_beat > 300:          # heartbeat for the Mac watcher's stall check
        last_beat = time.time()
        log(f"heartbeat: running {[r[0] for r in running]}, to do {len(todo)}")
    while todo and len(running) < args.par:
        j = todo.pop(0)
        jid = os.path.basename(j)[:-5]
        wd = os.path.join(TMP, jid)
        cmd = ["python3", os.path.join(HERE, "dna_check.py"), j, RES, "--threads", str(tpj), "--work", wd]
        if args.max_pairs:
            cmd += ["--max-pairs", str(args.max_pairs)]
        out = open(os.path.join(TMP, jid + ".out"), "a")
        running.append((jid, wd, subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT), time.time()))
        log("start", jid)
    for item in list(running):
        jid, wd, p, t0 = item
        if p.poll() is None:
            continue
        running.remove(item)
        res = os.path.join(RES, jid + ".json")
        if p.returncode == 0 and os.path.exists(res):
            r = json.load(open(res))
            v = {k: r[k]["verdict"] for k in r if isinstance(r[k], dict) and "verdict" in r[k]}
            log("done", jid, f"{time.time()-t0:.0f}s", json.dumps(v)[:300])
            shutil.rmtree(wd, ignore_errors=True)
        else:
            log("FAILED", jid, "rc", p.returncode, "(see tmp_dna/%s.out)" % jid)
    time.sleep(5)
log("dna driver finished")
