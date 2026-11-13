"""Job 1, step 3: orthologue-length consistency, fusion signal and presence in other
assemblies, by DIAMOND search of 1,000-aa windows of every giant.

Databases (data/refdb, not committed):
  ref  = UniProtKB/Swiss-Prot (all lengths) + UniProtKB unreviewed entries >= 3,000 aa
         (release 2026_04); used for orthologues and for the fusion test
  self = the 552 unique giant sequences themselves; used for near-identical copies in
         other assemblies
Per query sequence (all measured):
  windows, windows_with_ref_hit
  best other-species homologue: accession, organism, length, windows hit, mean identity,
  length ratio (query / homologue)
  best same-species homologue in UniProtKB (other entry), length and coverage
  fusion test: greedy number of distinct reference proteins needed to cover 80 % of the
  hit windows, and whether >= 2 of them are each covered >= 70 % over their own length
  on disjoint query segments ("segments match different complete proteins")
  near-identical copies (>= 95 % identity over >= 90 % of the query) in the giant set
  from another contig-set prefix (other assembly)
Usage: 12_orthology.py [threads]
"""
import csv
import gzip
import os
import re
import subprocess
import sys
import time
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
R = os.path.join(D, "refdb")
T = sys.argv[1] if len(sys.argv) > 1 else "8"
W = 1000


def read_fasta(path):
    name, buf = None, []
    op = gzip.open if path.endswith(".gz") else open
    for line in op(path, "rt"):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


stage0 = {r["md5"]: r for r in csv.DictReader(open(os.path.join(D, "stage0_table.tsv")), delimiter="\t")}
seqs = dict(read_fasta(os.path.join(D, "fasta", "unique_ge35992.fasta")))
cpath = os.path.join(D, "fasta", "controls.fasta")
controls = dict(read_fasta(cpath)) if os.path.exists(cpath) else {}
# one query per locus (longest live sequence); isoforms of the same gene inherit its row.
# (A first run over all 510 sequences in DIAMOND default mode was stopped after 2.8 CPU-h.)
rep_of, reps = {}, {}
for m, r in sorted(stage0.items(), key=lambda kv: -int(kv[1]["length"])):
    if r["triage"] == "T-excl":
        continue
    k = r["locus_key"]
    reps.setdefault(k, m)
    rep_of[m] = reps[k]
query_ids = sorted(set(reps.values())) + list(controls)
for c in controls:
    rep_of[c] = c
allseq = dict(seqs, **controls)

# ---- databases -------------------------------------------------------------------
t0 = time.time()
ref_db = os.path.join(R, "ref.dmnd")
if not os.path.exists(ref_db):
    with open(os.path.join(R, "ref.fasta"), "w") as f:
        for p in ("uniprot_sprot.fasta.gz", "uniprotkb_ge3000.fasta.gz"):
            for line in gzip.open(os.path.join(R, p), "rt"):
                f.write(line)
    subprocess.run(["diamond", "makedb", "--in", os.path.join(R, "ref.fasta"), "-d", ref_db[:-5],
                    "--threads", T, "--quiet"], check=True)
    os.remove(os.path.join(R, "ref.fasta"))
self_db = os.path.join(R, "self.dmnd")
if not os.path.exists(self_db):
    with open(os.path.join(R, "self.fasta"), "w") as f:
        for m, s in allseq.items():
            f.write(f">{m}\n{s}\n")
    subprocess.run(["diamond", "makedb", "--in", os.path.join(R, "self.fasta"), "-d", self_db[:-5], "--quiet"],
                   check=True)
qpath = os.path.join(R, "windows.fasta")
with open(qpath, "w") as f:
    for m in query_ids:
        s = allseq[m]
        for off in range(0, len(s), W):
            piece = s[off:off + W]
            if len(piece) >= 200:
                f.write(f">{m}|{off}\n{piece}\n")
print(f"databases ready {time.time() - t0:.0f}s", flush=True)

fmt = ["6", "qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue",
       "bitscore", "slen", "stitle"]
for name, db, extra in (("ref", ref_db, ["--max-target-seqs", "25", "--fast", "--max-hsps", "1"]),
                        ("self", self_db, ["--max-target-seqs", "50", "--fast"])):
    out = os.path.join(R, f"hits_{name}.tsv")
    if not os.path.exists(out):
        t1 = time.time()
        subprocess.run(["diamond", "blastp", "-q", qpath, "-d", db, "-o", out, "--threads", T, "--evalue", "1e-10",
                        "--outfmt", *fmt, "--quiet", *extra], check=True)
        print(f"diamond {name}: {time.time() - t1:.0f}s wall", flush=True)


def organism(title):
    m = re.search(r" OS=(.*?) OX=", title)
    return m.group(1) if m else ""


# ---- aggregate --------------------------------------------------------------------
nwin = {m: len(range(0, len(allseq[m]), W)) for m in query_ids}
ref = defaultdict(lambda: defaultdict(list))      # query -> target -> [(win, pident, qs, qe, ss, se)]
tinfo = {}
for line in open(os.path.join(R, "hits_ref.tsv")):
    f = line.rstrip("\n").split("\t")
    m, off = f[0].split("|")
    off = int(off)
    t = f[1]
    ref[m][t].append((off, float(f[2]), off + int(f[4]), off + int(f[5]), int(f[6]), int(f[7])))
    tinfo[t] = (int(f[10]), organism(f[11]), f[11].split(" OS=")[0].split(" ", 1)[-1][:60])
selfh = defaultdict(lambda: defaultdict(list))
for line in open(os.path.join(R, "hits_self.tsv")):
    f = line.rstrip("\n").split("\t")
    m, off = f[0].split("|")
    if f[1] == m:
        continue
    selfh[m][f[1]].append((int(off), float(f[2]), int(f[3])))

rows = []
for m in query_ids:
    s0 = stage0.get(m, {})
    org = s0.get("organism", "") if s0 else "Homo sapiens"
    qlen = len(allseq[m])
    hits = ref.get(m, {})
    win_any = {h[0] for t in hits.values() for h in t}
    row = {"md5": m, "windows": nwin[m], "windows_with_ref_hit": len(win_any)}
    # per-target summaries
    summ = []
    for t, hs in hits.items():
        wins = {h[0] for h in hs}
        slen, torg, tname = tinfo[t]
        tcov = bytearray(slen + 1)
        for h in hs:
            a, b = sorted((h[4], h[5]))
            tcov[a:b + 1] = b"\x01" * (b - a + 1)
        same_species = torg.split(" (")[0].lower() == org.split(" (")[0].lower() or \
            (org and torg and torg.split()[:2] == org.split()[:2])
        identical = slen == qlen and all(h[1] >= 99.9 for h in hs) and len(wins) >= 0.9 * nwin[m]
        summ.append({"t": t, "wins": wins, "n": len(wins), "slen": slen, "org": torg, "name": tname,
                     "pid": sum(h[1] for h in hs) / len(hs), "tcov": sum(tcov) / slen,
                     "qspan": (min(h[2] for h in hs), max(h[3] for h in hs)), "same": same_species,
                     "identical": identical})
    other = sorted([x for x in summ if not x["same"]], key=lambda x: (-x["n"], -x["pid"]))
    same = sorted([x for x in summ if x["same"] and not x["identical"]], key=lambda x: (-x["n"], -x["pid"]))
    ident = [x for x in summ if x["identical"]]
    if other:
        b = other[0]
        row.update({"best_other_species_hit": b["t"].split("|")[1] if "|" in b["t"] else b["t"],
                    "best_other_species_org": b["org"], "best_other_species_name": b["name"],
                    "best_other_species_len": b["slen"], "best_other_species_windows": b["n"],
                    "best_other_species_pid": round(b["pid"], 1),
                    "length_ratio_query_over_best": round(qlen / b["slen"], 2),
                    "best_other_species_target_cov": round(b["tcov"], 2)})
        # orthologue consistency among other-species hits covering >= 50 % of the query
        cons = [x for x in other if x["n"] >= 0.5 * nwin[m]]
        row["n_other_species_hits_cov50"] = len(cons)
        row["n_other_species_hits_cov50_len_within_20pct"] = sum(1 for x in cons if 0.8 <= qlen / x["slen"] <= 1.25)
    if same:
        b = same[0]
        row.update({"best_same_species_other_entry": b["t"].split("|")[1] if "|" in b["t"] else b["t"],
                    "best_same_species_len": b["slen"], "best_same_species_windows": b["n"],
                    "best_same_species_pid": round(b["pid"], 1)})
    row["identical_entries_in_uniprotkb"] = len(ident)
    row["identical_entry_accessions"] = ";".join(x["t"].split("|")[1] if "|" in x["t"] else x["t"] for x in ident)[:200]
    # fusion test: greedy cover of hit windows by distinct targets; complete-protein segments
    remaining, chosen = set(win_any), []
    pool = sorted(summ, key=lambda x: -x["n"])
    while remaining and len(chosen) < 10:
        best = max(pool, key=lambda x: len(x["wins"] & remaining), default=None)
        if not best or not (best["wins"] & remaining):
            break
        chosen.append(best)
        remaining -= best["wins"]
        if len(remaining) <= 0.2 * len(win_any):
            break
    row["targets_to_cover_80pct_hit_windows"] = len(chosen)
    complete = [x for x in chosen if x["tcov"] >= 0.7 and max(1000, 0.2 * qlen) <= x["slen"] <= 0.6 * qlen]
    disjoint = []
    for x in sorted(complete, key=lambda x: x["qspan"]):
        if not disjoint or x["qspan"][0] > disjoint[-1]["qspan"][1] - 200:
            disjoint.append(x)
    row["complete_proteins_on_disjoint_segments"] = len(disjoint)
    row["disjoint_segments_cover_frac"] = round(sum(x["qspan"][1] - x["qspan"][0] + 1 for x in disjoint) / qlen, 3)
    row["fusion_segments"] = ";".join(f"{x['qspan'][0]}-{x['qspan'][1]}:{x['name'][:30]}({x['slen']}aa,{x['org'][:25]})"
                                      for x in disjoint)[:400]
    # near-identical copies in the giant set from other assemblies
    copies = []
    for t, hs in selfh.get(m, {}).items():
        wins = {h[0] for h in hs}
        if len(wins) >= 0.9 * nwin[m] and sum(h[1] for h in hs) / len(hs) >= 95 and \
                abs(len(allseq[t]) - qlen) <= 0.05 * qlen:
            copies.append(t)
    row["near_identical_copies_in_giant_set"] = len(copies)
    row["near_identical_copy_accessions"] = ";".join(stage0[c]["rep_accession"] for c in copies if c in stage0)[:300]
    rows.append(row)

# isoforms / other sequences of the same locus inherit the representative's row
byq = {r["md5"]: r for r in rows}
for m, rep in rep_of.items():
    if m not in byq and rep in byq:
        rows.append(dict(byq[rep], md5=m, orthology_from_locus_representative=stage0[rep]["rep_accession"]))
cols = []
for r in rows:
    for k in r:
        if k not in cols:
            cols.append(k)
with open(os.path.join(D, "orthology_checks.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
print("rows", len(rows))
