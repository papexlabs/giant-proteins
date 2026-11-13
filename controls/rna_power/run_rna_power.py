"""Power test for the RNA-seq 'contradicted' rule (it fired for no model in the audit).

Wrong gene models are built for well-covered real models and scored against the SAME reads
with the verdict logic of cloud/summarize.py:
  shift    a fraction f of the introns is moved by 15-60 bp (both ends, same offset, so the
           intron length is kept; offsets above the 10-bp window of the placement guard)
  skip     a fraction f of the internal exons is removed from the model (each removal merges
           two introns into one longer intron)
A wrong intron counts as supported only if the reads contain that exact junction: the
per-run results store unique counts for the annotated introns and for the 50 most frequent
other junctions of each run (>= 2 reads), so a wrong intron that matches a rare real junction
is missed (this can only make a wrong model look worse, i.e. it favours detection).
Testability of a wrong intron uses the depth of the original flanking exons (best run).
Models: the human titin control and the eukaryotic models with verdict 'full' and >= 100
testable introns (data/rnaseq_support.tsv). Fractions 0.05-1.0, 20 random replicates each.
Output: controls/rna_power/rna_power.tsv and summary.json
"""
import csv
import json
import os
import random
import statistics
import tarfile
from collections import Counter, defaultdict

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
OUT = os.path.dirname(os.path.abspath(__file__))
MR, FULL, PART, FLANK, MINLOC = 3, 0.90, 0.50, 10.0, 50
FRACS = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
REPS = 20
rng = random.Random(20261113)

sup = {r["rep_accession"]: r for r in csv.DictReader(open(os.path.join(ROOT, "data", "rnaseq_support.tsv")), delimiter="\t")}
chosen = ["NP_001254479.2"] + [a for a, r in sorted(sup.items(), key=lambda x: -int(x[1]["testable_introns"] or 0))
                               if r["verdict"] == "full" and int(r["testable_introns"] or 0) >= 100 and a != "NP_001254479.2"]
md5 = {sup[a]["md5"]: a for a in chosen}

# per-run data of the chosen models (both archives; a _nosat re-run replaces its original)
runs = defaultdict(dict)      # acc -> run -> (locus dict, all junctions of that run)
for arch in ("rnaseq_run_results.tar.gz", "followup_run_results.tar.gz"):
    p = os.path.join(ROOT, "data", arch)
    with tarfile.open(p) as t:
        for m in t.getmembers():
            base = os.path.basename(m.name)
            if not m.name.endswith(".json") or base.startswith("._") or "/_" in m.name or base.endswith(".skip.json"):
                continue
            d = json.load(t.extractfile(m))
            loci = [L for L in d.get("loci", []) if L["md5"] in md5]
            if not loci:
                continue
            junc = Counter()
            for X in d["loci"]:
                for i in X["introns"]:
                    junc[(i[0], i[1])] = max(junc[(i[0], i[1])], i[2])
                for nj in X.get("novel_junctions", []):
                    junc[(nj[0], nj[1])] = max(junc[(nj[0], nj[1])], nj[2])
            nosat = "_nosat" in m.name
            for L in loci:
                acc = md5[L["md5"]]
                if d["run"] in runs[acc] and not nosat:
                    continue
                runs[acc][d["run"]] = (L, junc)


def verdict(introns, flank, uniq, leave, order="v2"):
    """cloud/summarize.py logic; introns = list of (U, M) summed over runs. order v1 = the rules
    used for the audit until 13 Nov; v2 = the 'under-supported where covered' test first."""
    n = len(introns)
    supported = [u >= MR or (u + m) >= MR for u, m in introns]
    counted = sum(supported)
    testable = [k for k in range(n) if flank[k] >= FLANK]
    t_uns = sum(1 for k in testable if not supported[k])
    ann = sum(u + m for u, m in introns)
    if uniq < MINLOC:
        return "no-expression"
    if (len(testable) >= 10 and t_uns >= 0.5 * len(testable)) or (leave > ann and leave >= 10):
        return "contradicted"
    if order == "v2" and len(testable) >= 10 and t_uns >= 0.2 * len(testable):
        return "weak, under-supported where covered"
    if counted >= FULL * n:
        return "full"
    if counted >= PART * n:
        return "partial"
    if len(testable) >= 10 and t_uns >= 0.2 * len(testable):
        return "weak, under-supported where covered"
    return "weak, reads too few"


def score(acc, model_introns, flank_of):
    """model_introns: list of (start, end); flank_of: function intron index -> best flank depth."""
    counts = []
    for k, (a, b) in enumerate(model_introns):
        u = m = 0
        for L, junc in runs[acc].values():
            ann = {(i[0], i[1]): (i[2], i[3]) for i in L["introns"]}
            if (a, b) in ann:
                u += ann[(a, b)][0]
                m += ann[(a, b)][1]
            else:
                u += junc.get((a, b), 0)
        counts.append((u, m))
    uniq = sum(L["reads_unique"] for L, _ in runs[acc].values())
    leave = sum(x[2] for L, _ in runs[acc].values() for x in L.get("novel_junctions", []) if x[3])
    fl_ = [flank_of(k) for k in range(len(model_introns))]
    return {o: verdict(counts, fl_, uniq, leave, o) for o in ("v1", "v2")}


rows = []
for acc in chosen:
    if not runs[acc]:
        continue
    L0 = next(iter(runs[acc].values()))[0]
    introns = [(i[0], i[1]) for i in L0["introns"]]
    exons = sorted((e[0], e[1]) for e in L0["exons"])
    depth = defaultdict(float)
    for L, _ in runs[acc].values():
        for e in L["exons"]:
            depth[(e[0], e[1])] = max(depth[(e[0], e[1])], e[3])
    left = {e[1]: depth[e] for e in exons}
    right = {e[0]: depth[e] for e in exons}
    fl = [min(left.get(a - 1, 0.0), right.get(b + 1, 0.0)) for a, b in introns]
    true_v = score(acc, introns, lambda k: fl[k])["v2"]
    n = len(introns)
    for kind in ("shift", "skip"):
        for f in FRACS:
            out = {"v1": Counter(), "v2": Counter()}
            for _ in range(REPS):
                if kind == "shift":
                    pick = set(rng.sample(range(n), max(1, round(f * n))))
                    model, flank = [], []
                    for k, (a, b) in enumerate(introns):
                        if k in pick:
                            s = rng.choice([-1, 1]) * rng.randint(15, 60)
                            model.append((a + s, b + s))
                        else:
                            model.append((a, b))
                        flank.append(fl[k])
                else:
                    internal = list(range(1, n))          # exon k sits between intron k-1 and k
                    pick = set(rng.sample(internal, max(1, round(f * len(internal)))))
                    model, flank = [], []
                    k = 0
                    cur_start, cur_left = introns[0][0], fl[0]
                    while k < n:
                        a, b = introns[k]
                        if k == 0 or (k not in pick):
                            cur_start = a
                        end = b
                        while k + 1 < n and (k + 1) in pick:
                            k += 1
                            end = introns[k][1]
                        model.append((cur_start, end))
                        # flank depth: exon before the merged intron, exon after it
                        lf = left.get(cur_start - 1, 0.0)
                        rf = right.get(end + 1, 0.0)
                        flank.append(min(lf, rf))
                        k += 1
                v = score(acc, model, lambda i: flank[i])
                for o in ("v1", "v2"):
                    out[o][v[o]] += 1
            rows.append({"model": acc, "organism": sup[acc]["organism"], "introns": n,
                         "testable_true": sum(1 for x in fl if x >= FLANK), "true_verdict": true_v, "error": kind,
                         "fraction_wrong": f, "replicates": REPS,
                         **{f"{o}_{k2}": out[o].get(k2, 0) for o in ("v1", "v2")
                            for k2 in ("contradicted", "full", "partial", "weak, under-supported where covered",
                                       "weak, reads too few")}})
    print(acc, sup[acc]["organism"], n, true_v, flush=True)

cols = list(rows[0].keys())
with open(os.path.join(OUT, "rna_power.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
summ = {"models": len({r["model"] for r in rows}), "replicates": REPS, "fractions": FRACS, "by_kind_fraction": {}}
for kind in ("shift", "skip"):
    for f in FRACS:
        rr = [r for r in rows if r["error"] == kind and r["fraction_wrong"] == f]
        tot = sum(r["replicates"] for r in rr)
        summ["by_kind_fraction"][f"{kind} {f}"] = {o: {
            "contradicted": round(sum(r[f"{o}_contradicted"] for r in rr) / tot, 3),
            "under-supported where covered": round(sum(r[f"{o}_weak, under-supported where covered"] for r in rr) / tot, 3),
            "partial": round(sum(r[f"{o}_partial"] for r in rr) / tot, 3),
            "full": round(sum(r[f"{o}_full"] for r in rr) / tot, 3)} for o in ("v1", "v2")}
summ["true_models_verdicts"] = dict(Counter(r["true_verdict"] for r in rows if r["error"] == "shift" and r["fraction_wrong"] == FRACS[0]))
json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
