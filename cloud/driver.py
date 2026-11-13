#!/usr/bin/env python3
"""RNA-seq junction-support driver for the giant-protein audit, v2 (runs on the box or,
for testing, on the Mac). One scheduler keeps three kinds of work going at once:

  stage A  (per run, up to --par runs in parallel, all mates of a run in parallel)
           stream ENA FASTQ over HTTPS -> gzip -dc -> head (read cap) ->
           minimap2 -x sr --sam-hit-only against the CDS EXONS (+-100 bp) of the
           species' giant loci (bundle/exons/<species>.fa; no genome needed) ->
           sam2fq.py (candidate cap). Only reads that touch an exon are kept.
  prep     (per species, background) download the genome assembly, place the loci on
           it (prep_genome.py) and build the whole-genome minimap2 splice:sr index.
           The whole genome is the decoy: a candidate read counts for a giant only if
           its best (primary) spliced alignment in the genome lands on the giant.
  stage B  (per run, as soon as its candidates and its species' index are ready)
           minimap2 -x splice:sr -t N -> SAM file -> samtools view -F 0x904 -L loci.bed
           (only alignments that touch a locus) -> rnaseq_locus.py -> results/<sp>/<run>.json

Durability: every finished run is written at once (atomic rename) to results/, which the
Mac syncs every 10 minutes; runs with a result are skipped after a relaunch. Runs that
cannot be processed get results/<sp>/<run>.skip.json with the reason. driver.log is
appended to (launch.sh pushes the synced log back before starting).

Per-run limits (logged as stop_rule in each result):
  saturation  stop when the candidate reads reach --sat-cov x (exon bp / 100) per mate
  read cap    at most --max-reads reads per FASTQ file (first reads of the file)
  time cap    the download of each FASTQ file stops after --run-time-limit seconds
              (curl --max-time); everything read so far is aligned and used
Deadline: no new run starts within 40 min of the deadline (/etc/gp-deadline on the box).

--download (added 13 Nov 2026 after ENA streams broke with curl rc=18 in earlier runs): each
FASTQ file is first downloaded to disk with resumed retries (curl -C -, up to 12 attempts, a
stalled transfer is cut after 120 s below 100 kB/s and resumed), checked against ENA's byte
size (fastq_sizes in the bundle), and only then decompressed and filtered from the file. The
file is deleted with the run's other files once the run's result has been written. The genome
is downloaded with up to 6 attempts and checked with gzip -t. A run starts only if the disk
keeps >= 30 GB free after its files (counting the full size of downloads still in progress).

Usage: driver.py <bundle_dir> <work_dir> [--threads N] [--par K] [--max-reads R]
                 [--sat-cov C] [--run-time-limit S] [--deadline-epoch T] [--no-shutdown]
                 [--max-species N] [--max-runs N] [--download] [--only SPECIES_KEY]
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import time

UA = "PapexResearchBot/1.0 (+https://papex.org)"
HERE = os.path.dirname(os.path.abspath(__file__))

ap = argparse.ArgumentParser()
ap.add_argument("bundle")
ap.add_argument("work")
ap.add_argument("--threads", type=int, default=os.cpu_count())
ap.add_argument("--par", type=int, default=8, help="runs in stage A at the same time (2 streams each)")
ap.add_argument("--max-reads", type=int, default=60_000_000, help="reads per FASTQ file")
ap.add_argument("--sat-cov", type=float, default=60.0, help="candidate cap = sat-cov x exon_bp/100 per mate")
ap.add_argument("--run-time-limit", type=int, default=1800)
ap.add_argument("--deadline-epoch", type=float, default=0)
ap.add_argument("--max-species", type=int, default=10 ** 9)
ap.add_argument("--max-runs", type=int, default=10 ** 9, help="per species (testing)")
ap.add_argument("--no-shutdown", action="store_true")
ap.add_argument("--download", action="store_true", help="download each FASTQ file to disk (resumed retries), then map")
ap.add_argument("--only", default="", help="process only this species key (testing)")
args = ap.parse_args()

RES = os.path.join(args.work, "results")
TMP = os.path.join(args.work, "tmp")
os.makedirs(RES, exist_ok=True)
os.makedirs(TMP, exist_ok=True)
LOG = open(os.path.join(args.work, "driver.log"), "a")


def log(*a):
    msg = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + " ".join(str(x) for x in a)
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()


def mem_bytes():
    try:
        return int(next(l for l in open("/proc/meminfo") if l.startswith("MemTotal")).split()[1]) * 1024
    except (OSError, StopIteration):
        out = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True).stdout.strip()
        return int(out or 0)


def write_json(path, obj):
    with open(path + ".tmp", "w") as f:
        json.dump(obj, f, separators=(",", ":"))
    os.replace(path + ".tmp", path)


deadline = args.deadline_epoch
if not deadline and os.path.exists("/etc/gp-deadline"):
    deadline = float(open("/etc/gp-deadline").read().strip())
MEM = mem_bytes()
# minimap2 splice:sr index: building peaks at ~16 bytes per genome base (measured on the
# Mac: 242-Mb chr2 -> 3.85 GB RSS); mapping with the saved index needs ~5-6 bytes per base.
B_PAR = 2 if MEM >= 100e9 else 1                  # concurrent stage-B jobs (each loads an index)
PREP_PAR = 2 if MEM >= 200e9 else 1               # concurrent genome preparations (index builds)
B_THREADS = max(4, args.threads // 4)
A_THREADS = max(1, (args.threads - B_PAR * B_THREADS) // (2 * args.par))
log(f"driver v2 start: threads {args.threads}, par {args.par}, mem {MEM/1e9:.0f} GB, stage-B jobs {B_PAR} "
    f"x {B_THREADS} threads, stage-A minimap2 threads per mate {A_THREADS}, max reads/file {args.max_reads}, "
    f"sat-cov {args.sat_cov}, run time limit {args.run_time_limit}s, deadline "
    f"{time.strftime('%FT%TZ', time.gmtime(deadline)) if deadline else 'none'}")
for p in [os.path.join(dp, f) for dp, _, fs in os.walk(RES) for f in fs if f.endswith((".tmp", ".part"))]:
    os.remove(p)                                   # partial results of a killed box

# ---------------------------------------------------------------------------------------
jobs = []
with open(os.path.join(args.bundle, "jobs.tsv")) as f:
    cols = f.readline().rstrip("\n").split("\t")
    for line in f:
        jobs.append(dict(zip(cols, line.rstrip("\n").split("\t"))))
if args.only:
    jobs = [j for j in jobs if j["species_key"] == args.only]
jobs = jobs[:args.max_species]
species = {}
queue = []                                         # global run queue in job order
for j in jobs:
    sp = j["species_key"]
    os.makedirs(os.path.join(RES, sp), exist_ok=True)
    if os.path.exists(os.path.join(RES, sp, "DONE")):
        continue
    gbp = int(j.get("genome_bp") or 0)
    if MEM and gbp * 16 > 0.8 * MEM:
        log("SKIP species (genome too large for this box's memory):", sp, gbp)
        write_json(os.path.join(RES, sp, "DONE"), {"skipped": "genome too large for memory"})
        continue
    runs = json.loads(j["runs_json"])[:args.max_runs]
    exon_fa = os.path.join(args.bundle, "exons", f"{sp}.fa")
    exon_bp = sum(len(l.strip()) for l in open(exon_fa) if not l.startswith(">")) if os.path.exists(exon_fa) else 0
    species[sp] = {"job": j, "runs": runs, "exon_fa": exon_fa, "exon_bp": exon_bp, "prep": None,
                   "prep_state": "todo", "wd": os.path.join(TMP, sp)}
    for r in runs:
        if os.path.exists(os.path.join(RES, sp, r["run"] + ".json")) or \
                os.path.exists(os.path.join(RES, sp, r["run"] + ".skip.json")):
            continue
        queue.append((sp, r))
log(f"{len(species)} species to do, {len(queue)} runs in the queue")


def species_pending(sp):
    return [r for r in species[sp]["runs"]
            if not os.path.exists(os.path.join(RES, sp, r["run"] + ".json"))
            and not os.path.exists(os.path.join(RES, sp, r["run"] + ".skip.json"))]


def start_prep(sp):
    S = species[sp]
    os.makedirs(S["wd"], exist_ok=True)
    wd = S["wd"]
    if args.download:
        get = (f"for t in 1 2 3 4 5 6; do curl -sfL --retry 8 --retry-delay 5 --connect-timeout 30 -A '{UA}' "
               f"-o {wd}/genome.fna.gz '{S['job']['genome_url']}' && gzip -t {wd}/genome.fna.gz && break; "
               f"rm -f {wd}/genome.fna.gz; sleep $((t * 30)); done; test -s {wd}/genome.fna.gz; ")
    else:
        get = (f"curl -sfL --retry 8 --retry-delay 5 --connect-timeout 30 -A '{UA}' -o {wd}/genome.fna.gz "
               f"'{S['job']['genome_url']}'; ")
    cmd = (f"set -e; {get}"
           f"python3 {HERE}/prep_genome.py {wd}/genome.fna.gz {args.bundle}/loci/{sp}.json {wd} > {wd}/prep.json; "
           f"minimap2 -x splice:sr -t 8 -I 64G -d {wd}/genome.mmi {wd}/genome.fna.gz 2> {wd}/index.log; "
           f"rm -f {wd}/genome.fna.gz; touch {wd}/READY")
    S["prep"] = subprocess.Popen(cmd, shell=True, executable="/bin/bash")
    S["prep_state"] = "running"
    S["prep_t0"] = time.time()
    log("prep start", sp, S["job"]["genome_url"].split("/")[-1], f"{int(S['job'].get('genome_bp') or 0)/1e9:.2f} Gb")


def start_stage_a(sp, r):
    S = species[sp]
    wd = S["wd"]
    os.makedirs(wd, exist_ok=True)
    if not os.path.exists(os.path.join(wd, "exons.mmi")):
        subprocess.run(f"minimap2 -x sr -d {wd}/exons.mmi {S['exon_fa']} 2> /dev/null", shell=True, check=True,
                       executable="/bin/bash")
    urls = [u for u in r.get("fastq_urls", []) if u]
    cap = int(args.sat_cov * S["exon_bp"] / 100 / max(1, len(urls))) if S["exon_bp"] else 0
    parts = []
    sizes = r.get("fastq_sizes") or []
    for i, url in enumerate(urls, 1):
        base = f"{wd}/{r['run']}_{i}"
        if args.download:
            exp = int(sizes[i - 1]) if i - 1 < len(sizes) and sizes[i - 1] else 0
            f = f"{base}.fastq.gz"
            # rc file: 0 = complete file on disk (size checked when ENA gives it), else curl's last rc
            parts.append(
                f"( ok=0; rc=1; for t in 1 2 3 4 5 6 7 8 9 10 11 12; do "
                f"curl -sfL -C - --connect-timeout 30 --speed-limit 100000 --speed-time 120 -A '{UA}' -o {f} '{url}'; rc=$?; "
                f"sz=$(stat -c %s {f} 2>/dev/null || stat -f %z {f} 2>/dev/null || echo 0); "
                f"if [ $rc -eq 0 ] && {{ [ {exp} -eq 0 ] || [ $sz -eq {exp} ]; }}; then ok=1; break; fi; "
                f"if [ $rc -eq 33 ] || {{ [ {exp} -gt 0 ] && [ $sz -gt {exp} ]; }}; then rm -f {f}; fi; "
                f"echo \"attempt $t rc=$rc size=$sz of {exp}\" >> {base}.dl.log; sleep $((t * 15)); done; "
                f"if [ $ok -eq 1 ]; then "
                f"gzip -dc {f} 2>/dev/null | head -n {4 * args.max_reads} "
                f"| minimap2 -t {A_THREADS} -a -x sr --secondary=no --sam-hit-only {wd}/exons.mmi - 2> {base}.mm.log "
                f"| python3 {HERE}/sam2fq.py {i} {cap} > {base}.cand.fq 2> {base}.cap.log; "
                f"g=${{PIPESTATUS[0]}}; if [ $g -eq 0 ] || [ $g -eq 141 ]; then echo 0 > {base}.rc; else echo gz$g > {base}.rc; fi; "
                f"else echo $rc > {base}.rc; fi ) &")
            continue
        # the time cap is applied to the download only (curl --max-time): the stream then ends
        # normally and every stage downstream flushes what it has (no buffered reads lost).
        # No curl --retry: a retried stream would start again from byte 0 and duplicate reads;
        # a run that yields no reads at all is re-queued instead (at most 2 more attempts).
        parts.append(
            f"( curl -sfL --connect-timeout 30 --max-time {args.run_time_limit} "
            f"-A '{UA}' '{url}' | gzip -dc 2>/dev/null | head -n {4 * args.max_reads} "
            f"| minimap2 -t {A_THREADS} -a -x sr --secondary=no --sam-hit-only {wd}/exons.mmi - 2> {base}.mm.log "
            f"| python3 {HERE}/sam2fq.py {i} {cap} > {base}.cand.fq 2> {base}.cap.log; "
            f"echo ${{PIPESTATUS[0]}} > {base}.rc ) &")
    cmd = " ".join(parts) + " wait"
    return {"sp": sp, "r": r, "proc": subprocess.Popen(cmd, shell=True, executable="/bin/bash"),
            "t0": time.time(), "n_mates": len(urls), "cap": cap}


def stage_a_stats(a):
    wd = species[a["sp"]]["wd"]
    reads, cands, rules = 0, 0, []
    for i in range(1, a["n_mates"] + 1):
        base = f"{wd}/{a['r']['run']}_{i}"
        mate_reads = 0
        try:
            mate_reads = sum(int(m) for m in re.findall(r"mapped (\d+) sequences", open(base + ".mm.log").read()))
        except OSError:
            pass
        reads += mate_reads
        try:
            cands += sum(1 for _ in open(base + ".cand.fq")) // 4
        except OSError:
            pass
        rc = open(base + ".rc").read().strip() if os.path.exists(base + ".rc") else "?"
        capped = os.path.exists(base + ".cap.log") and "cap of" in open(base + ".cap.log").read()
        # rc = curl's exit status: 0 complete, 28 time cap (--max-time), 23/141 stopped by a
        # downstream cap (head / sam2fq closed the pipe)
        if capped:
            rules.append("saturation")
        elif mate_reads >= args.max_reads * 0.999:
            rules.append("read cap")
        elif rc == "28":
            rules.append("time cap")
        elif rc.startswith("gz"):
            rules.append(f"gzip error rc={rc[2:]}")
        elif rc not in ("0", "23", "141"):
            rules.append(f"stream error curl rc={rc}")
        else:
            rules.append("complete")
    return reads, cands, rules


def start_stage_b(sp, r, a_info):
    S = species[sp]
    wd = S["wd"]
    run = r["run"]
    out = os.path.join(RES, sp, run + ".json")
    cmd = (f"set -o pipefail; cat {wd}/{run}_*.cand.fq | minimap2 -t {B_THREADS} -a -x splice:sr {wd}/genome.mmi - "
           f"2> {wd}/{run}.B.log > {wd}/{run}.sam && "
           f"samtools view -c -F 0x904 {wd}/{run}.sam > {wd}/{run}.total && "
           f"samtools view -F 0x904 -L {wd}/loci.bed {wd}/{run}.sam "
           f"| python3 {HERE}/rnaseq_locus.py {wd}/loci_found.json {run} {wd}/{run}.total > {out}.part")
    return {"sp": sp, "r": r, "a": a_info, "out": out,
            "proc": subprocess.Popen(cmd, shell=True, executable="/bin/bash"), "t0": time.time()}


def finish_run(b):
    sp, r, a = b["sp"], b["r"], b["a"]
    wd = species[sp]["wd"]
    part = b["out"] + ".part"
    if b["proc"].returncode != 0 or not os.path.exists(part) or os.path.getsize(part) == 0:
        log("  run FAILED stage B", r["run"], "rc", b["proc"].returncode)
        write_json(os.path.join(RES, sp, r["run"] + ".skip.json"), {"run": r["run"], "reason": "stage B failed"})
    else:
        res = json.load(open(part))
        res.update({"candidate_reads": a["cands"], "reads_processed": a["reads"], "stop_rule": a["rules"],
                    "stage_a_s": round(a["secs"]), "stage_b_s": round(time.time() - b["t0"]),
                    "run_meta": {k: r[k] for k in r if k not in ("fastq_urls", "_try")}})
        write_json(b["out"], res)
        os.remove(part)
        log("  run", r["run"], f"A {a['secs']:.0f}s B {time.time()-b['t0']:.0f}s reads {a['reads']} cand {a['cands']}",
            "/".join(sorted(set(a["rules"]))), "; ".join(
                f"{L['rep_accession']}:{L['reads_unique']}u {L['introns_with_ge3_exact_reads']}"
                f"(+multi {L['introns_with_ge3_incl_multimapped']})/{L['n_introns']}i" for L in res["loci"])[:300])
    for f in os.listdir(wd):
        if f.startswith(r["run"] + "_") or f.startswith(r["run"] + "."):
            os.remove(os.path.join(wd, f))


active_a, ready_b, active_b = [], [], []
last_beat = 0
stop_new = False
while True:
    now = time.time()
    if deadline and now > deadline - 40 * 60 and not stop_new:
        stop_new = True
        log("STOP RULE: within 40 min of the deadline, no new runs or species are started")
    # skip runs without FASTQ URLs (logged, durable)
    for item in [q for q in queue if not [u for u in q[1].get("fastq_urls", []) if u]]:
        queue.remove(item)
        sp, r = item
        write_json(os.path.join(RES, sp, r["run"] + ".skip.json"), {"run": r["run"], "reason": "no FASTQ URL in ENA"})
        log("  run SKIPPED (no FASTQ URL)", sp, r["run"])
    # genome preparation: for the species of the next runs in line (lookahead)
    upcoming = []
    for sp, _ in [(a["sp"], 0) for a in active_a] + [(x["sp"], 0) for x in ready_b] + queue[:args.par * 2]:
        if sp not in upcoming:
            upcoming.append(sp)
    for sp in upcoming:
        S = species[sp]
        if S["prep_state"] == "todo" and not stop_new and \
                sum(1 for x in species.values() if x["prep_state"] == "running") < PREP_PAR:
            start_prep(sp)
    for sp, S in species.items():
        if S["prep_state"] == "running" and S["prep"].poll() is not None:
            ok = S["prep"].returncode == 0 and os.path.exists(os.path.join(S["wd"], "READY"))
            S["prep_state"] = "ready" if ok else "failed"
            info = open(os.path.join(S["wd"], "prep.json")).read().strip() if os.path.exists(
                os.path.join(S["wd"], "prep.json")) else ""
            log("prep", "ready" if ok else "FAILED", sp, f"{time.time()-S['prep_t0']:.0f}s", info[:200])
            if not ok:
                for r in species_pending(sp):
                    write_json(os.path.join(RES, sp, r["run"] + ".skip.json"),
                               {"run": r["run"], "reason": "genome preparation failed"})
    # stage A: fill free slots in queue order (species whose preparation failed are skipped)
    while queue and len(active_a) < args.par and not stop_new:
        sp, r = queue[0]
        if species[sp]["prep_state"] == "failed":
            queue.pop(0)
            continue
        if args.download:
            need = sum(int(x or 0) for x in (r.get("fastq_sizes") or [])) or int(r.get("bytes") or 0)
            busy = sum(sum(int(x or 0) for x in (a["r"].get("fastq_sizes") or [])) or int(a["r"].get("bytes") or 0)
                       for a in active_a)
            if shutil.disk_usage(args.work).free - busy < need + 30e9:
                if not active_a and not active_b and not ready_b:
                    queue.pop(0)
                    write_json(os.path.join(RES, sp, r["run"] + ".skip.json"), {"run": r["run"], "reason": "not enough disk"})
                    log("  run SKIPPED (not enough disk)", sp, r["run"], f"{need/1e9:.1f} GB")
                    continue
                break                                  # wait for running work to free disk
        queue.pop(0)
        active_a.append(start_stage_a(sp, r))
    for a in list(active_a):
        if a["proc"].poll() is not None:
            active_a.remove(a)
            reads, cands, rules = stage_a_stats(a)
            if reads == 0 and any(x.startswith(("stream error", "gzip error")) for x in rules) and a["r"].get("_try", 0) < 2:
                a["r"]["_try"] = a["r"].get("_try", 0) + 1
                queue.append((a["sp"], a["r"]))
                log("  run stream failed, re-queued", a["r"]["run"], rules, "attempt", a["r"]["_try"] + 1)
                for f in os.listdir(species[a["sp"]]["wd"]):
                    if f.startswith(a["r"]["run"] + "_"):
                        os.remove(os.path.join(species[a["sp"]]["wd"], f))
                continue
            ready_b.append({"sp": a["sp"], "r": a["r"],
                            "a": {"reads": reads, "cands": cands, "rules": rules, "secs": time.time() - a["t0"]}})
    # stage B: runs whose species index is ready
    for x in list(ready_b):
        if len(active_b) >= B_PAR:
            break
        st = species[x["sp"]]["prep_state"]
        if st == "ready":
            ready_b.remove(x)
            active_b.append(start_stage_b(x["sp"], x["r"], x["a"]))
        elif st == "failed":
            ready_b.remove(x)
    for b in list(active_b):
        if b["proc"].poll() is not None:
            active_b.remove(b)
            finish_run(b)
    # species completion: all runs have a result or a skip record
    for sp, S in species.items():
        if S.get("done"):
            continue
        busy = any(x["sp"] == sp for x in active_a + ready_b + active_b) or any(q[0] == sp for q in queue)
        if not busy and (not species_pending(sp) or S["prep_state"] == "failed" or stop_new):
            if not species_pending(sp):
                write_json(os.path.join(RES, sp, "DONE"), {"finished": time.strftime("%FT%TZ", time.gmtime())})
                log("species done", sp)
            if S["prep"] and S["prep"].poll() is None:
                S["prep"].kill()
            shutil.rmtree(S["wd"], ignore_errors=True)
            S["done"] = True
    if now - last_beat > 600:
        last_beat = now
        la = os.getloadavg()[0] if hasattr(os, "getloadavg") else 0
        log(f"heartbeat: load {la:.1f}, stage A {len(active_a)}, waiting for index {len(ready_b)}, "
            f"stage B {len(active_b)}, preps {sum(1 for x in species.values() if x['prep_state']=='running')}, "
            f"queue {len(queue)}")
    if not queue and not active_a and not ready_b and not active_b and \
            not any(S["prep_state"] == "running" for S in species.values()):
        break
    if stop_new and not active_a and not active_b and not ready_b:
        break
    time.sleep(2)

log("driver finished")
if os.path.exists("/etc/gp-deadline") and not args.no_shutdown:
    open(os.path.join(args.work, "ALL_DONE"), "w").write(time.strftime("%FT%TZ", time.gmtime()))
    log("powering off (instance-initiated shutdown = terminate)")
    subprocess.Popen("sleep 600; sudo shutdown -h now", shell=True)
