#!/usr/bin/env python3
"""Build the RNA-seq job bundle (cloud/bundle/) from the Job 1 locus checks.

Selection (eukaryotes only; bacteria/MAGs have no usable RNA-seq):
  - positive controls: PKZILLA-1 (KAL1508129.1) and human titin (curated RefSeq)
  - every eukaryotic sequence with provisional grade B or C (RNA-seq can move it to
    A/B or D); grade A loci from RefSeq Gnomon evidence and grade D loci are skipped,
    except one negative control (the 98,182-aa sheep-hybrid repeat model)
Per species: genome assembly URL (NCBI FTP over HTTPS, via elink nuccore->assembly),
up to K RNA-seq runs from ENA (Illumina, RNA-Seq, 0.2-6 GB, different studies first;
muscle/heart samples first when the locus is titin-like).
Order: controls, then species whose best locus is > 45,212 aa, then by locus length,
so the run can stop at any point with the most informative species done.
Usage: make_jobs.py [K=4] [--select repeat-flag] [--out DIR] [--suffix S]
  --select repeat-flag   only the models never sent to the RNA stage because they were graded
                         D1 under the first repeat rule (now C with a repeat flag; rnaseq_verdict
                         'not-tested (graded D1 ...') -> cloud/bundle_rep_rna, species keys + '_rep'
  --out, --suffix        bundle directory and species-key suffix (defaults: cloud/bundle, '')
Each run also records its per-file FASTQ sizes (fastq_sizes), which the driver's --download mode
uses to check that a file on disk is complete.
"""
import csv
import gzip
import json
import os
import re
import sys
import urllib.parse

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from http_util import get  # noqa: E402

D = os.path.join(ROOT, "data")
_a = sys.argv[1:]


def _opt(name, default):
    if name in _a:
        i = _a.index(name)
        v = _a[i + 1]
        del _a[i:i + 2]
        return v
    return default


SELECT = _opt("--select", "")
OUT = os.path.abspath(_opt("--out", os.path.join(ROOT, "cloud", "bundle_rep_rna" if SELECT == "repeat-flag" else "bundle")))
SUFFIX = _opt("--suffix", "_rep" if SELECT == "repeat-flag" else "")
os.makedirs(os.path.join(OUT, "loci"), exist_ok=True)
K = int(_a[0]) if _a else 4      # runs per species; top species K + 2
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
SINGLE_CELL = re.compile(r"scRNA|single[- ]?cell|snRNA|10x|drop-?seq", re.I)
MUSCLE = re.compile(r"muscle|heart|cardiac|myo|skeletal|ventric|atri", re.I)

checks = list(csv.DictReader(open(os.path.join(D, "locus_checks.tsv")), delimiter="\t"))
specs = {}
for p in ("specs.jsonl", "control_specs.jsonl"):
    for line in open(os.path.join(D, "loci", p)):
        s = json.loads(line)
        specs[s["md5"]] = s
avail = {r["organism"]: r for r in csv.DictReader(open(os.path.join(D, "species_data_availability.tsv")), delimiter="\t")}
avail["Homo sapiens"] = {"taxid": "9606"}

sel = []
for r in checks:
    s = specs.get(r["md5"])
    if not s or not s["segments"] or r["superkingdom"] != "Eukaryota":
        continue
    if SELECT == "repeat-flag":
        if r["rnaseq_verdict"].startswith("not-tested (graded D1"):
            sel.append((r, s, "repeat-flag"))
        continue
    ctrl = r["rep_accession"] in ("KAL1508129.1",) or r["md5"].startswith("control_")
    neg = r["rep_accession"] == "KAI4551101.1"
    fusion_check = r["grade"] == "D" and "D5" in r["grade_reason"]     # merged-model check
    if ctrl or neg or fusion_check or r["grade"] in ("B", "C"):
        sel.append((r, s, "positive-control" if ctrl else "negative-control" if neg else
                    "fusion-check" if fusion_check else "audit"))

CACHE = os.path.join(OUT, "assembly_cache.json")
acache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}


def assembly_url(contig):
    if contig in acache:
        return acache[contig]
    url, size = "", 0
    try:
        r = json.loads(get(EU + "esearch.fcgi", {"db": "nuccore", "term": f"{contig}[accn]", "retmode": "json"}))
        ids = r["esearchresult"]["idlist"]
        if ids:
            l = json.loads(get(EU + "elink.fcgi", {"dbfrom": "nuccore", "db": "assembly", "id": ids[0],
                                                   "retmode": "json"}))
            links = [x for ls in l["linksets"] for db in ls.get("linksetdbs", []) for x in db["links"]]
            if links:
                s = json.loads(get(EU + "esummary.fcgi", {"db": "assembly", "id": links[0], "retmode": "json"}))
                d = s["result"][s["result"]["uids"][0]]
                ftp = d["ftppath_refseq"] if contig.startswith(("NC_", "NW_", "NZ_")) and d.get("ftppath_refseq") \
                    else d["ftppath_genbank"]
                if ftp:
                    url = ftp.replace("ftp://", "https://") + "/" + ftp.rstrip("/").split("/")[-1] + "_genomic.fna.gz"
                    size = int(d.get("meta", "").split('total_length" ')[-1].split("<")[0].split(">")[-1] or 0) \
                        if 'total_length' in d.get("meta", "") else 0
    except Exception as e:  # noqa: BLE001
        print("assembly lookup failed", contig, e)
    if not url and SELECT:
        url, size = ena_assembly(contig)
    acache[contig] = (url, size)
    json.dump(acache, open(CACHE, "w"))
    return acache[contig]


def ena_assembly(contig):
    """Fallback for INSDC records NCBI has not linked to an assembly yet (recent ENA
    submissions, e.g. OZ chromosomes): the contig's ENA project -> the assembly of that
    project whose record lists the contig (chromosome accession or WGS prefix). Genome
    FASTA from the ENA browser API (gzip); size estimated from the largest ENA assembly
    of the same taxon with a base count (the new assemblies have none yet)."""
    try:
        t = get(f"https://www.ebi.ac.uk/ena/browser/api/embl/{contig}", {"lineLimit": 12}).decode()
        pr = re.search(r"^PR\s+Project:(PRJ\w+)", t, re.M)
        tx = re.search(r'/db_xref="taxon:(\d+)"', get(f"https://www.ebi.ac.uk/ena/browser/api/embl/{contig}",
                                                         {"lineLimit": 80}).decode())
        if not pr:
            return "", 0
        rows = get("https://www.ebi.ac.uk/ena/portal/api/search",
                   {"result": "assembly", "query": f'study_accession="{pr.group(1)}"',
                    "fields": "accession,version", "format": "tsv"}).decode().splitlines()[1:]
        for line in rows:
            acc, ver = line.split("\t")[:2]
            x = get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{acc}.{ver}").decode()
            wgs = re.search(r"<PREFIX>(\w+)</PREFIX>", x)
            if f'accession="{contig}"' in x or (wgs and contig.startswith(wgs.group(1))):
                size = 0
                if tx:
                    bc = get("https://www.ebi.ac.uk/ena/portal/api/search",
                             {"result": "assembly", "query": f"tax_eq({tx.group(1)})", "fields": "base_count",
                              "format": "tsv"}).decode().splitlines()[1:]
                    size = max([int(v) for b in bc for v in b.split("\t") if v.isdigit()] + [0])
                return (f"https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}.{ver}?download=true&gzip=true", size)
    except Exception as e:  # noqa: BLE001
        print("ENA assembly lookup failed", contig, e)
    return "", 0


# Controls: the PKZILLA-1 genome strain (12B1) has only two RNA-seq runs (48-58 GB); one is
# added on top of the size-capped picks so the positive control can work. Human: skeletal muscle.
FORCE_RUNS = {"97485": ["SRR23517916"]}
TISSUE_FILTER = {"9606": ' AND tissue_type="skeletal muscle"'}


def ena_runs(taxid, titin_like, k):
    q = f'tax_eq({taxid}) AND library_strategy="RNA-Seq" AND instrument_platform="ILLUMINA"'
    q += TISSUE_FILTER.get(taxid, "")
    fields = "run_accession,fastq_ftp,fastq_bytes,read_count,library_layout,study_accession,sample_title,tissue_type"
    t = get("https://www.ebi.ac.uk/ena/portal/api/search",
            {"result": "read_run", "query": q, "fields": fields, "format": "tsv", "limit": 2000}).decode()
    rows = list(csv.DictReader(t.splitlines(), delimiter="\t"))
    good, forced = [], []
    for x in rows:
        if not x.get("fastq_ftp"):
            continue
        b = sum(int(v) for v in x["fastq_bytes"].split(";") if v)
        if x["run_accession"] in FORCE_RUNS.get(taxid, []):
            forced.append((x, b))
            continue
        if not 0.2e9 <= b <= 6e9 or not fastq_urls(x):
            continue
        if SINGLE_CELL.search(x.get("sample_title", "")):   # 3'-tag single-cell: poor junction coverage
            continue
        good.append((x, b))
    if not good and SELECT:                             # only larger runs exist: allow up to 15 GB
        for x in rows:                                  # (downloaded to disk, not streamed)
            b = sum(int(v) for v in x["fastq_bytes"].split(";") if v) if x.get("fastq_ftp") else 0
            if 6e9 < b <= 15e9 and fastq_urls(x) and not SINGLE_CELL.search(x.get("sample_title", "")):
                good.append((x, b))
    good.sort(key=lambda xb: (barcode_like(xb[0]),
                              -(titin_like and bool(MUSCLE.search(xb[0]["sample_title"] + " " +
                                                                  xb[0].get("tissue_type", "")))),
                              -int(xb[0]["read_count"] or 0)))
    picked, studies = list(forced), {}
    for x, b in good:                                   # round-robin over studies
        studies.setdefault(x["study_accession"], []).append((x, b))
    while len(picked) < k + len(forced) and any(studies.values()):
        for st in list(studies):
            if studies[st] and len(picked) < k + len(forced):
                picked.append(studies[st].pop(0))
    return [{"run": x["run_accession"], "bytes": b, "read_count": int(x["read_count"] or 0),
             "layout": x["library_layout"], "study": x["study_accession"],
             "sample": (x["sample_title"] + " " + x.get("tissue_type", ""))[:80],
             "fastq_urls": fastq_urls(x), "fastq_sizes": fastq_sizes(x)}
            for x, b in picked], len(rows), len(good)


def fastq_sizes(x):
    """Byte sizes of the files fastq_urls() returns, in the same order (0 when ENA gives none)."""
    ftp = x["fastq_ftp"].split(";")
    sizes = [int(v) for v in x["fastq_bytes"].split(";")] if x.get("fastq_bytes") else []
    by = {"https://" + u: (sizes[i] if i < len(sizes) else 0) for i, u in enumerate(ftp)}
    return [by.get(u, 0) for u in fastq_urls(x)]


def fastq_urls(x):
    """Mate files _1/_2 when ENA has them; otherwise every FASTQ file listed (some runs
    registered as PAIRED have a single file); every file is read as single-end.
    Single-cell style runs (one mate much smaller: a barcode/UMI read of ~28 bp) keep only
    the cDNA mate."""
    ftp = x["fastq_ftp"].split(";")
    sizes = [int(v) for v in x["fastq_bytes"].split(";")] if x.get("fastq_bytes") else []
    files = [(u, sizes[i] if i < len(sizes) else 0) for i, u in enumerate(ftp) if u.endswith(".fastq.gz")]
    mates = [f for f in files if f[0].endswith(("_1.fastq.gz", "_2.fastq.gz"))] or files
    if len(mates) == 2 and min(m[1] for m in mates) > 0 and max(m[1] for m in mates) > 1.8 * min(m[1] for m in mates):
        mates = [max(mates, key=lambda m: m[1])]
    return ["https://" + u for u, _ in mates]


def barcode_like(x):
    sizes = [int(v) for v in x["fastq_bytes"].split(";") if v] if x.get("fastq_bytes") else []
    return len(sizes) == 2 and min(sizes) > 0 and max(sizes) > 1.8 * min(sizes)


IDX = json.load(open(os.path.join(D, "loci", "dna_index.json")))
_dna = {}


def region_seq(s):
    for lid, v in IDX.items():
        if s["md5"] in v["md5s"] and v["how"] not in ("fail", "fail-too-large"):
            if lid not in _dna:
                _dna[lid] = "".join(l.strip() for l in gzip.open(os.path.join(D, "loci", "dna", lid + ".fa.gz"), "rt")
                                    if not l.startswith(">"))
            return _dna[lid], v["start"]
    return "", 0


def exon_records(s, pad=100):
    """CDS exons +-pad bp (merged) as + strand sequences: the stage-A filter index."""
    seq, start = region_seq(s)
    if not seq:
        return []
    iv = sorted((max(1, a - pad), b + pad) for a, b in s["segments"])
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return [(f"{s['md5'][:12]}_{a}", seq[max(0, a - start):b - start + 1]) for a, b in merged]


def probe_for(s):
    """300 bp of + strand genomic sequence at the first CDS base; lets the driver find the
    locus when the assembly FASTA names the contig differently (e.g. CM vs WGS alias)."""
    lo = min(a for a, b in s["segments"])
    for lid, v in IDX.items():
        if s["md5"] in v["md5s"] and v["how"] not in ("fail", "fail-too-large"):
            if lid not in _dna:
                _dna[lid] = "".join(l.strip() for l in gzip.open(os.path.join(D, "loci", "dna", lid + ".fa.gz"), "rt")
                                    if not l.startswith(">"))
            return _dna[lid][lo - v["start"]:lo - v["start"] + 300], lo
    return "", lo


groups = {}
for r, s, role in sel:
    # one group per organism and assembly: groups are merged later if two contig prefixes
    # resolve to the same assembly (the key collision of the first run lost 3 loci)
    groups.setdefault((r["organism"], s["contig"].split(".")[0][:6]), []).append((r, s, role))
_by_key = {}
jobs = []
for (org, _), items in groups.items():
    taxid = avail.get(org, {}).get("taxid", "")
    if not taxid:
        continue
    url, gsize = assembly_url(items[0][1]["contig"])
    if not url:
        print("no assembly URL for", org)
        continue
    titin_like = any("Ig/FN3" in r["panel_class"] for r, _, _ in items)
    roles0 = {role for _, _, role in items}
    top = bool(roles0 & {"positive-control", "fusion-check"}) or max(int(r["length"]) for r, _, _ in items) > 45212
    runs, n_all, n_ok = ena_runs(taxid, titin_like, K + 2 if top else K)
    key = re.sub(r"[^A-Za-z0-9]+", "_", org)[:40] + "_" + url.split("/")[-1].split("?")[0].split("_")[1][:12] + SUFFIX
    if key in _by_key:                     # same assembly as an earlier group: merge the loci
        prev = _by_key[key]
        items = prev["items"] + items
        jobs = [j for j in jobs if j["species_key"] != key]
    _by_key[key] = {"items": items}
    loci = []
    for r, s, _ in items:
        pr, pp = probe_for(s)
        loci.append({"locus_id": s["md5"][:12], "md5": s["md5"], "rep_accession": s["rep_accession"],
                     "contig": s["contig"], "strand": s["strand"], "segments": s["segments"],
                     "grade0": r["grade"], "probe": pr, "probe_pos": pp})
    json.dump(loci, open(os.path.join(OUT, "loci", f"{key}.json"), "w"))
    os.makedirs(os.path.join(OUT, "exons"), exist_ok=True)
    with open(os.path.join(OUT, "exons", f"{key}.fa"), "w") as fe:
        for _, s, _ in items:
            for name, sq in exon_records(s):
                fe.write(f">{name}\n{sq}\n")
    roles = {role for _, _, role in items}
    maxlen = max(int(r["length"]) for r, _, _ in items)
    gb = sum(x["bytes"] for x in runs) / 1e9
    est = 5 + (gsize / 1e9) * 6 + gb / 0.08 / 60          # min: index ~6 min/Gb, stream ~80 MB/s
    grades = {r["grade"] for r, _, _ in items}
    prio = (0 if "positive-control" in roles else
            1 if "repeat-flag" in roles and maxlen > 45212 else
            1 if maxlen > 45212 or "fusion-check" in roles else
            2 if "negative-control" in roles else
            3 if "C" in grades else 4, -maxlen)
    jobs.append({"species_key": key, "organism": org, "taxid": taxid, "role": ",".join(sorted(roles)),
                 "n_loci": len(loci), "max_locus_aa": maxlen, "genome_url": url, "genome_bp": gsize,
                 "ena_rnaseq_runs_total": n_all, "ena_runs_0.2-6GB": n_ok, "n_runs": len(runs),
                 "runs_gb": round(gb, 2), "est_wall_min": round(est, 1), "prio": prio,
                 "runs_json": json.dumps(runs, separators=(",", ":"))})
    print(f"{org[:30]:30s} loci {len(loci):3d} runs {len(runs)}/{n_ok}/{n_all} {gb:6.1f} GB genome {gsize/1e9:.2f} Gb")
jobs = [j for j in jobs if j["n_runs"] > 0]
jobs.sort(key=lambda j: j["prio"])
cols = [c for c in jobs[0] if c != "prio"]
with open(os.path.join(OUT, "jobs.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for j in jobs:
        f.write("\t".join(str(j[c]) for c in cols) + "\n")
empty = [(j["species_key"], r["run"]) for j in jobs for r in json.loads(j["runs_json"]) if not r["fastq_urls"]]
print("runs without FASTQ URL:", len(empty), empty[:5])
print("jobs", len(jobs), "runs", sum(j["n_runs"] for j in jobs), "GB", round(sum(j["runs_gb"] for j in jobs), 1),
      "est wall h (sequential)", round(sum(j["est_wall_min"] for j in jobs) / 60, 1))
