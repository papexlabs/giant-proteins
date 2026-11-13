#!/usr/bin/env python3
"""Pre-launch check of cloud/bundle: every species has loci, an exon FASTA and a genome
URL; every run has at least one HTTPS FASTQ URL. With --head, every FASTQ URL and genome
URL is checked with an HTTP HEAD request (generic User-Agent, no credentials) and the
real byte sizes are summed. Prints the data volume the driver will read under its caps.
Exit code 1 if anything is missing.
Usage: check_bundle.py [--head] [--max-reads 60000000] [--bundle bundle_followup_rna]
"""
import csv
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HERE, sys.argv[sys.argv.index("--bundle") + 1]) if "--bundle" in sys.argv else os.path.join(HERE, "bundle")
UA = "PapexResearchBot/1.0 (+https://papex.org)"
head = "--head" in sys.argv
max_reads = int(sys.argv[sys.argv.index("--max-reads") + 1]) if "--max-reads" in sys.argv else 60_000_000

jobs = list(csv.DictReader(open(os.path.join(B, "jobs.tsv")), delimiter="\t"))
problems, urls = [], []
total_bytes = capped_bytes = 0
for j in jobs:
    sp = j["species_key"]
    for p in (os.path.join(B, "loci", sp + ".json"), os.path.join(B, "exons", sp + ".fa")):
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            problems.append(f"{sp}: missing {os.path.basename(p)}")
    if not j["genome_url"].startswith("https://"):
        problems.append(f"{sp}: genome URL {j['genome_url']!r}")
    urls.append(j["genome_url"])
    for r in json.loads(j["runs_json"]):
        u = [x for x in r.get("fastq_urls", []) if x]
        if not u or not all(x.startswith("https://") for x in u):
            problems.append(f"{sp} {r['run']}: no usable FASTQ URL")
        urls += u
        total_bytes += r["bytes"]
        # read cap per file: estimate bytes per read from the run's own size / read count
        per_file_reads = r["read_count"] or 1
        bpr = r["bytes"] / max(1, per_file_reads * max(1, len(u)))
        capped_bytes += min(r["bytes"], int(max_reads * bpr * max(1, len(u))))


def head_size(u):
    try:
        req = urllib.request.Request(u, method="HEAD", headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=60) as resp:
            return u, resp.status, int(resp.headers.get("Content-Length") or 0)
    except Exception as e:  # noqa: BLE001
        return u, repr(e)[:80], 0


if head:
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(head_size, urls))
    bad = [r for r in res if r[1] != 200]
    for u, st, _ in bad:
        problems.append(f"HEAD {st}: {u}")
    print(f"HEAD checked {len(res)} URLs, {len(bad)} failed, {sum(r[2] for r in res)/1e9:.1f} GB in total")
print(f"{len(jobs)} species, {sum(len(json.loads(j['runs_json'])) for j in jobs)} runs; FASTQ volume "
      f"{total_bytes/1e9:.0f} GB, after the {max_reads/1e6:.0f}M-reads-per-file cap about {capped_bytes/1e9:.0f} GB")
for p in problems:
    print("PROBLEM", p)
sys.exit(1 if problems else 0)
