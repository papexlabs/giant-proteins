"""Per-species availability of independent evidence data, for the live giant entries.

For each organism (taxid) carrying at least one live, non-synthetic giant sequence:
  - SRA experiments: RNA-seq (short read), long-read RNA (PacBio/ONT, transcriptomic),
    long-read genomic (PacBio/ONT, WGS strategy)
  - median run size (GB) and bases of up to 20 sampled RNA-seq runs (esummary)
  - PRIDE Archive projects matching the organism name (keyword search, rough)
No API keys, no contact parameters; NCBI <= 2.5 req/s.
Output: data/species_data_availability.tsv
"""
import csv
import json
import os
import re
import statistics
import sys
import urllib.parse

sys.path.insert(0, os.path.dirname(__file__))
from http_util import get  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

uniq = list(csv.DictReader(open(os.path.join(D, "giant_unique_sequences.tsv")), delimiter="\t"))
ncbi = {r["accession"]: r for r in csv.DictReader(open(os.path.join(D, "ncbi_protein_ge35992.tsv")), delimiter="\t")}
upc = {r["upi"]: r for r in csv.DictReader(open(os.path.join(D, "uniparc_ge35992.tsv")), delimiter="\t")}

taxa = {}
for r in uniq:
    if r["currently_live"] != "yes" or r["patent_or_synthetic"] == "yes":
        continue
    tid, org = "", r["organism"]
    for a in r["ncbi_accessions"].split(";"):
        if a in ncbi:
            tid = ncbi[a]["taxid"]
            break
    if not tid and r["upi"]:
        tid = upc[r["upi"].split(";")[0]]["taxid"]
    if not tid:
        continue
    t = taxa.setdefault(tid, {"taxid": tid, "organism": org, "n_giant_seqs": 0, "max_len": 0,
                              "superkingdom": r["superkingdom"]})
    t["n_giant_seqs"] += 1
    t["max_len"] = max(t["max_len"], int(r["length"]))


def count(term):
    r = json.loads(get(EU + "esearch.fcgi", {"db": "sra", "term": term, "retmode": "json", "retmax": 20},
                       timeout=30))
    return int(r["esearchresult"]["count"]), r["esearchresult"]["idlist"]


CACHE = os.path.join(D, "raw", "species_availability_cache.json")
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
rows = []
for tid, t in sorted(taxa.items(), key=lambda kv: -kv[1]["max_len"]):
    if tid in cache:
        t.update(cache[tid])
        rows.append(t)
        continue
    org = f"txid{tid}[Organism:noexp]"
    n_rna, ids = count(f'{org} AND "rna seq"[Strategy]')
    n_lr_rna, _ = count(f'{org} AND ("pacbio smrt"[Platform] OR "oxford nanopore"[Platform]) '
                        f'AND (transcriptomic[Source] OR "rna seq"[Strategy] OR "FL-cDNA"[Strategy])')
    n_lr_dna, _ = count(f'{org} AND ("pacbio smrt"[Platform] OR "oxford nanopore"[Platform]) AND genomic[Source]')
    sizes, bases = [], []
    if ids:
        s = json.loads(get(EU + "esummary.fcgi", {"db": "sra", "id": ",".join(ids), "retmode": "json"},
                           timeout=30))
        for uid in s["result"].get("uids", []):
            runs = s["result"][uid].get("runs", "")
            for m in re.finditer(r'total_spots="(\d+)" total_bases="(\d+)"', runs):
                bases.append(int(m.group(2)))
            expxml = s["result"][uid].get("expxml", "")
            m = re.search(r'total_size="(\d+)"', expxml + runs)
            if m:
                sizes.append(int(m.group(1)))
    pride = ""
    try:
        name = t["organism"].split(" x ")[0]
        q = urllib.parse.quote(name)
        p = json.loads(get(f"https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects?keyword={q}&pageSize=100",
                           timeout=30, retries=2))
        pride = sum(1 for x in p if any(name.lower() in (o or "").lower()
                                        for o in (x.get("organisms") or []) + [x.get("title", "")]))
    except Exception as e:  # noqa: BLE001
        pride = "error"
    t.update({
        "sra_rnaseq_experiments": n_rna, "sra_longread_rna_experiments": n_lr_rna,
        "sra_longread_genomic_experiments": n_lr_dna,
        "median_rnaseq_run_gbases": round(statistics.median(bases) / 1e9, 2) if bases else "",
        "median_rnaseq_run_size_gb": round(statistics.median(sizes) / 1e9, 2) if sizes else "",
        "pride_projects_name_match": pride,
    })
    rows.append(t)
    cache[tid] = {k: t[k] for k in ("sra_rnaseq_experiments", "sra_longread_rna_experiments",
                                     "sra_longread_genomic_experiments", "median_rnaseq_run_gbases",
                                     "median_rnaseq_run_size_gb", "pride_projects_name_match")}
    json.dump(cache, open(CACHE, "w"))
    print(tid, t["organism"], n_rna, n_lr_rna, n_lr_dna, t["median_rnaseq_run_size_gb"], pride, flush=True)

cols = list(rows[0].keys())
with open(os.path.join(D, "species_data_availability.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
