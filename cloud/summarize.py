#!/usr/bin/env python3
"""Aggregate the per-run RNA-seq results (v1 and v2 driver formats) into one row per giant
sequence model and assign a transcript-support verdict.

Per intron, reads are summed over all runs of the species:
  U = spliced reads matching the annotated intron exactly, unique (MAPQ >= 10)
  M = the same, multi-mapped (MAPQ < 10)
  supported (unique)      U >= MIN_READS
  supported (multi only)  U < MIN_READS and U + M >= MIN_READS
  none                    U + M < MIN_READS
Multi-only support counts towards the verdict only for models WITHOUT a repeat flag
(Stage 1 rules D1/D2/D3: translated tandem repeat, non-coding DNA, ORF open in other
frames): inside a tandem-repeat model a multi-mapped read cannot tell copies apart.
Coverage breadth = share of CDS-exon bases covered by >= 1 unique read in at least one
run (per exon the best run; a lower bound of the union).

Verdicts (thresholds set on the controls, see notes/RESULTS.md):
  no-rnaseq-runs the species had no usable Illumina RNA-seq run in ENA (0.2-6 GB, not single-cell)
  not-tested     no run result for this model (reason given)
  no-expression  fewer than MIN_LOCUS_READS unique reads over all runs
  full           counted-supported introns >= FULL_FRAC of all introns
  partial        >= PARTIAL_FRAC
  testable intron: both flanking exons reach mean depth >= FLANK_DEPTH (10) in at least one
                 run; only there does a missing junction count against the model (revised
                 12 Nov after the titin follow-up: coverage of these giants is often 3'-heavy,
                 so a gene-wide "expected reads per intron" over-called under-support)
  weak           "under-supported where covered" if >= 20 % of >= 10 testable introns lack
                 support (since 13 Nov checked before full/partial: the power test in
                 controls/rna_power showed that models with 10-40 % wrong introns were otherwise
                 called "partial"); else < PARTIAL_FRAC supported: "under-supported where covered" if >= 20 % of >= 10
                 testable introns lack support, else "unsupported introns lie where reads are
                 too few"
  contradicted   >= 10 testable introns and >= 50 % of them unsupported, or junctions that
                 leave the gene outnumber the reads on all annotated junctions (>= 10 reads)
  single-exon    no introns: coverage breadth only
Runs: one result per (model, run). A run re-done without the saturation stop (species key
suffix _nosat, follow-up 12 Nov) replaces the earlier result of the same run; loci re-run after
the job-key collision come from species keys with suffix _extra. Planned models and job lists
are read from every bundle given (default: cloud/bundle and cloud/bundle_followup_rna).
Placement guard: a model relocated between non-identical sequences is counted at lifted
coordinates (data/followups/relocation_liftover.json, scripts/21_relocation_liftover.py); any
other model whose unsupported introns show up in the reads 1-10 bp away at both ends (>= 10
such introns) is "not-tested (model coordinates shifted against the reads)".
Usage: summarize.py [--results cloud_out/results] [--out data/rnaseq_support.tsv]
"""
import argparse
import csv
import glob
import json
import os
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--results", default=os.path.join(ROOT, "cloud_out", "results"))
ap.add_argument("--bundle", action="append", default=None,
                help="bundle dir(s); default cloud/bundle and cloud/bundle_followup_rna")
ap.add_argument("--checks", default=os.path.join(ROOT, "data", "locus_checks.tsv"))
ap.add_argument("--out", default=os.path.join(ROOT, "data", "rnaseq_support.tsv"))
ap.add_argument("--min-reads", type=int, default=3, help="reads for an intron to count as supported")
ap.add_argument("--min-locus-reads", type=int, default=50, help="unique reads below which a model is 'no-expression'")
ap.add_argument("--full-frac", type=float, default=0.90)
ap.add_argument("--partial-frac", type=float, default=0.50)
ap.add_argument("--contra-frac", type=float, default=0.10)
ap.add_argument("--verdict-order", choices=("v1", "v2"), default="v2",
                help="v2: the 'under-supported where covered' test runs before full/partial (13 Nov)")
ap.add_argument("--flank-depth", type=float, default=10.0,
                help="mean depth both flanking exons must reach (in one run) for an intron to be testable")
args = ap.parse_args()
# bundle_rep_rna (repeat-flag models, 13 Nov) counts as planned only once results of that stage exist
BUNDLES = args.bundle or [os.path.join(ROOT, "cloud", b) for b in ("bundle", "bundle_followup_rna", "bundle_rep_rna")
                          if os.path.isdir(os.path.join(ROOT, "cloud", b)) and
                          (b != "bundle_rep_rna" or glob.glob(os.path.join(args.results, "*_rep", "*.json")))]
SUFFIXES = ("_nosat", "_extra", "_rep")
# repeat-flag models of the _rep stage (13 Nov): which introns have a junction sequence unique in
# the model (controls/rna_rep/junction_mappability.json, computed before any result existed).
# Pre-registered guard (cloud/bundle_rep_rna/RULES.txt): 'contradicted' moves such a model to D only
# if it also holds on the introns with a unique junction (or through junctions leaving the gene).
MAPPABILITY = {}
_mp = os.path.join(ROOT, "controls", "rna_rep", "junction_mappability.json")
if os.path.exists(_mp):
    MAPPABILITY = {x["md5"]: x["per_intron_genomic_order"] for x in json.load(open(_mp)) if x.get("md5")
                   and x.get("per_intron_genomic_order")}


def base_key(sp):
    for x in SUFFIXES:
        if sp.endswith(x):
            return sp[:-len(x)]
    return sp

checks = {r["md5"]: r for r in csv.DictReader(open(args.checks), delimiter="\t")} if os.path.exists(args.checks) else {}

# models sent to the RNA stage (bundle loci files) -> species key
planned = {}
for b in BUNDLES:
    for p in glob.glob(os.path.join(b, "loci", "*.json")):
        sp = base_key(os.path.basename(p)[:-5])
        for L in json.load(open(p)):
            planned[L["md5"]] = (sp, L["rep_accession"])

# models placed on another sequence whose intron coordinates moved under a proper lift-over
# (scripts/21_relocation_liftover.py): junction reads are counted at the lifted coordinates
LIFT = {}
lp = os.path.join(ROOT, "data", "followups", "relocation_liftover.json")
if os.path.exists(lp):
    LIFT = json.load(open(lp))
agg = defaultdict(lambda: {"runs": 0, "runs_with_reads": 0, "uniq": 0, "multi": 0, "reads_processed": 0,
                           "introns": None, "exons": None, "leave": 0, "alt": 0, "species": "", "rep": "",
                           "v1_runs": 0})
# one result per (model, run): a _nosat re-run replaces the saturation-stopped original
best = {}
for p in sorted(glob.glob(os.path.join(args.results, "*", "*.json"))):
    if p.endswith(".skip.json") or os.path.basename(os.path.dirname(p)).startswith("_"):
        continue
    d = json.load(open(p))
    sp_dir = os.path.basename(os.path.dirname(p))
    run = d.get("run") or os.path.basename(p)[:-5]
    prio = 0 if sp_dir.endswith("_nosat") else 1
    for L in d.get("loci", []):
        k = (L["md5"], run)
        if k not in best or prio < best[k][0]:
            best[k] = (prio, sp_dir, d, L)
superseded = len([1 for p in glob.glob(os.path.join(args.results, "*", "*.json"))
                  if not p.endswith(".skip.json") and not os.path.basename(os.path.dirname(p)).startswith("_")]) - \
    len({(v[1], k[1]) for k, v in best.items()})
for (md5, run), (prio, sp_dir, d, L) in sorted(best.items(), key=lambda x: (x[0][0], x[0][1])):
    sp = base_key(sp_dir)
    a = agg[L["md5"]]
    a["species"], a["rep"] = sp, L["rep_accession"]
    a.setdefault("nosat_runs", 0)
    a["nosat_runs"] += prio == 0
    a["runs"] += 1
    a["v1_runs"] += "stop_rule" not in d
    a["runs_with_reads"] += L["reads_unique"] > 0
    a["uniq"] += L["reads_unique"]
    a["multi"] += L["reads_multi"]
    a["reads_processed"] += d.get("reads_processed") or 0
    a["leave"] += sum(x[2] for x in L.get("novel_junctions", []) if x[3])
    if a["introns"] is None:
        a["introns"] = [[i[0], i[1], 0, 0, 0] for i in L["introns"]]
        a["exons"] = [[e[0], e[1], 0.0, 0.0] for e in L["exons"]]
    if L["md5"] in LIFT:
        # junction reads of this run at reference coordinates: annotated junctions of every
        # model in the run (unique reads, the best model per junction) and novel junctions
        per = {}
        for X in d.get("loci", []):
            for i in X["introns"]:
                per[(i[0], i[1])] = max(per.get((i[0], i[1]), 0), i[2])
            for nj in X.get("novel_junctions", []):
                per[(nj[0], nj[1])] = max(per.get((nj[0], nj[1]), 0), nj[2])
        for k, x in enumerate(LIFT[L["md5"]]["introns_lifted"]):
            a["introns"][k][2] += per.get(tuple(x), 0) if x else 0
    else:
        for k, i in enumerate(L["introns"]):
            a["introns"][k][2] += i[2]
            a["introns"][k][3] += i[3]
            a["introns"][k][4] += i[4] + i[5]
    nov = a.setdefault("novel", defaultdict(int))
    for nj in L.get("novel_junctions", []):
        nov[(nj[0], nj[1])] += nj[2]
    for k, e in enumerate(L["exons"]):
        a["exons"][k][2] = max(a["exons"][k][2], e[2])
        a["exons"][k][3] += e[3]
    # per intron: best (over runs) of the smaller depth of its two flanking exons
    left = {e[1]: e[3] for e in L["exons"]}
    right = {e[0]: e[3] for e in L["exons"]}
    fl = a.setdefault("flank", [0.0] * len(L["introns"]))
    for k, i in enumerate(L["introns"]):
        fl[k] = max(fl[k], min(left.get(i[0] - 1, 0.0), right.get(i[1] + 1, 0.0)))

job_species = set()
for b in BUNDLES:
    if os.path.exists(os.path.join(b, "jobs.tsv")):
        job_species |= {base_key(r["species_key"]) for r in csv.DictReader(open(os.path.join(b, "jobs.tsv")), delimiter="\t")}
jp = os.path.join(BUNDLES[0], "jobs.tsv")
# models the job builder selected (eukaryotic B/C, fusion checks, controls) but whose loci file
# was overwritten: two job groups of one species got the same species key (fixed in make_jobs)
eligible = {m for m, c in checks.items() if c.get("superkingdom") == "Eukaryota" and
            (c.get("grade") in ("B", "C") or "D5 " in c.get("grade_reason", "") or c.get("control") or
             c.get("rep_accession") == "KAI4551101.1")}
lost = {m for m in eligible if m not in planned and m not in agg}
SEGS = {}
sp_path = os.path.join(ROOT, "data", "loci", "specs.jsonl")
if os.path.exists(sp_path):
    for line in open(sp_path):
        x = json.loads(line)
        SEGS[x["md5"]] = len(x["segments"])
COLLIDED = set()
if os.path.exists(jp):
    keys = [(r["species_key"], r["organism"]) for r in csv.DictReader(open(jp), delimiter="\t")]
    COLLIDED = {o for k, o in keys if sum(1 for k2, _ in keys if k2 == k) > 1}

REPEAT_RULES = ("D1 ", "D2 ", "D3 ")
rows = []
for m in sorted(set(planned) | set(agg) | lost):
    a = agg.get(m)
    c = checks.get(m, {})
    sp, rep = planned.get(m, (a["species"] if a else "", a["rep"] if a else ""))
    repeat = any(x in c.get("grade_reason", "") for x in REPEAT_RULES)
    row = {"md5": m, "rep_accession": rep or c.get("rep_accession", ""), "organism": c.get("organism", ""),
           "length": c.get("length", ""), "species_key": sp, "stage1_grade": c.get("grade", ""),
           "repeat_flag": "yes" if repeat else "no"}
    if not a:
        if m in lost:
            if not SEGS.get(m):
                v = "not-tested (locus not resolved in Stage 1)"
            elif "repeat flag" in c.get("grade_reason", ""):
                v = "not-tested (graded D1 when the RNA stage was planned; D1 rule revised 13 Nov)"
            elif c.get("organism") in COLLIDED:
                v = "not-tested (job-list key collision)"
            else:
                v = "not-tested (no assembly link found)"
        elif sp and sp not in job_species:
            v = "no-rnaseq-runs (none usable in ENA)"
        else:
            v = "not-tested"
        row.update({"verdict": v})
        rows.append(row)
        continue
    MR = args.min_reads
    n = len(a["introns"])
    uq = sum(1 for i in a["introns"] if i[2] >= MR)
    mo = sum(1 for i in a["introns"] if i[2] < MR and i[2] + i[3] >= MR)
    none = n - uq - mo
    counted = uq + (0 if repeat else mo)
    ex_bp = sum(e[1] - e[0] + 1 for e in a["exons"])
    breadth = sum((e[1] - e[0] + 1) * e[2] for e in a["exons"]) / max(1, ex_bp)
    depth = sum((e[1] - e[0] + 1) * e[3] for e in a["exons"]) / max(1, ex_bp)
    annotated_reads = sum(i[2] + i[3] for i in a["introns"])
    # testable introns: both flanking exons reach depth >= FLANK in at least one run, so a
    # correct intron should show spliced reads there. Missing support only counts against a
    # model where it is testable (coverage of giant genes is often uneven, e.g. 3'-heavy).
    sup = [(i[2] >= MR) or (not repeat and i[2] + i[3] >= MR) for i in a["introns"]]
    testable = [k for k in range(n) if a["flank"][k] >= args.flank_depth]
    t_unsup = sum(1 for k in testable if not sup[k])
    # coordinate-shift signature: unsupported introns whose junction the reads show k bp away
    # at BOTH ends (1 <= |k| <= 10), i.e. the model is placed off the genome the reads map to
    shifted = 0
    if m not in LIFT:
        for i in a["introns"]:
            if i[2] < MR and any(a["novel"].get((i[0] + k, i[1] + k), 0) >= MR
                                 for k in range(-10, 11) if k):
                shifted += 1
    uj_check = {}
    if m in MAPPABILITY and n:
        per = MAPPABILITY[m]
        lens = [i[1] - i[0] + 1 for i in a["introns"]]
        flags = None
        for cand in (per, per[::-1]):
            if len(cand) == n and len({lens[k] - cand[k][0] for k in range(n)}) == 1:
                flags = [c[1] for c in cand]
                break
        if flags is None:
            uj_check = {"testable_unique_junction": "", "contradicted_on_unique_junctions": "not checked (intron lists differ)"}
        else:
            tu = [k for k in testable if flags[k]]
            tu_unsup = sum(1 for k in tu if not sup[k])
            c_u = (len(tu) >= 10 and tu_unsup >= 0.5 * len(tu)) or (a["leave"] > annotated_reads and a["leave"] >= 10)
            uj_check = {"introns_unique_junction": sum(flags), "testable_unique_junction": len(tu),
                        "testable_unique_junction_unsupported": tu_unsup,
                        "introns_with_ge3_reads_incl_multi": sum(1 for i in a["introns"] if i[2] + i[3] >= MR),
                        "contradicted_on_unique_junctions": "yes" if c_u else "no"}
    if a["uniq"] < args.min_locus_reads:
        verdict = "no-expression"
    elif shifted >= 10:
        verdict = f"not-tested (model coordinates shifted against the reads: {shifted} introns seen 1-10 bp away)"
    elif n == 0:
        verdict = "single-exon"
    elif (len(testable) >= 10 and t_unsup >= 0.5 * len(testable)) or (a["leave"] > annotated_reads and a["leave"] >= 10):
        verdict = "contradicted"
    elif args.verdict_order == "v2" and len(testable) >= 10 and t_unsup >= 0.2 * len(testable):
        # v2 (13 Nov, after the power test in controls/rna_power): a model with >= 20 % of its
        # covered introns unsupported is under-supported even if >= 50 % of all introns are
        # supported; v1 checked this only below 50 % and so called such models 'partial'
        verdict = "weak (under-supported where covered)"
    elif counted >= args.full_frac * n:
        verdict = "full"
    elif counted >= args.partial_frac * n:
        verdict = "partial"
    elif len(testable) >= 10 and t_unsup >= 0.2 * len(testable):
        verdict = "weak (under-supported where covered)"
    else:
        verdict = "weak (unsupported introns lie where reads are too few)"
    row.update({
        "runs": a["runs"], "runs_with_unique_reads": a["runs_with_reads"], "v1_runs": a["v1_runs"],
        "runs_without_saturation_stop": a.get("nosat_runs", 0),
        "reads_processed_v2": a["reads_processed"], "reads_unique": a["uniq"], "reads_multi": a["multi"],
        "n_introns": n, "introns_unique_support": uq, "introns_multi_only": mo, "introns_none": none,
        "introns_counted_supported": counted, "frac_introns_supported": round(counted / n, 3) if n else "",
        "cds_breadth": round(breadth, 3), "cds_mean_depth": round(depth, 1),
        "reads_on_annotated_junctions": annotated_reads, "reads_alt_partner_junctions": sum(i[4] for i in a["introns"]),
        "reads_junctions_leaving_gene": a["leave"], "testable_introns": len(testable),
        "testable_introns_unsupported": t_unsup,
        "introns_shifted_junction": shifted,
        "lifted_over": "yes: " + LIFT[m]["relocation"] if m in LIFT else "",
        "verdict": verdict,
    })
    row.update(uj_check)
    rows.append(row)

cols = []
for r in rows:
    for k in r:
        if k not in cols:
            cols.append(k)
os.makedirs(os.path.dirname(args.out), exist_ok=True)
with open(args.out, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: -int(r["length"] or 0)):
        f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
from collections import Counter  # noqa: E402
print(f"{len(rows)} models -> {args.out}; verdicts {dict(Counter(r['verdict'] for r in rows))}; "
      f"{len(best)} (model, run) results kept, {superseded} run file(s) superseded by a _nosat re-run")
