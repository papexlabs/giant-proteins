#!/usr/bin/env python3
"""Write DATA_MANIFEST.tsv: what the repository stores and what its scripts download instead.

Columns: kind (stored | downloaded | local-only), item (path or accession), bytes (stored: file size;
downloaded: size given by the archive, empty when the archive gives none), source, used_by.
- stored: every tracked file under data/, controls/, proteomics/, cloud/bundle*/ and paper/supplementary/
  (from `git ls-files`), with its size.
- downloaded: RNA-seq runs (cloud/bundle*/jobs.tsv, ENA fastq_bytes), genome assemblies (genome_url; size in bases),
  source reads of the DNA checks and controls (ENA accessions), proteomics files (PRIDE file lists).
- local-only: git-ignored folders that scripts fill (reads, genomes, raw spectra), listed so a rerun knows them.
Usage: python3 release/make_manifest.py   (no network)
"""
import csv
import glob
import json
import os
import subprocess

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
rows = []

tracked = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
for f in tracked:
    if f.startswith(("data/", "controls/", "proteomics/", "paper/supplementary/")) or (f.startswith("cloud/bundle")):
        if f.endswith((".py", ".sh")):
            continue
        rows.append(("stored", f, os.path.getsize(os.path.join(ROOT, f)), "this repository", ""))

seen = set()
for b in sorted(glob.glob(os.path.join(ROOT, "cloud", "bundle*", "jobs.tsv"))):
    stage = os.path.basename(os.path.dirname(b))
    for j in csv.DictReader(open(b), delimiter="\t"):
        g = j["genome_url"]
        if g and (g, stage) not in seen:
            seen.add((g, stage))
            rows.append(("downloaded", g, "", f"genome assembly, {int(j.get('genome_bp') or 0):,} bases (NCBI or ENA)",
                         f"cloud/driver.py ({stage}, {j['species_key']})"))
        for r in json.loads(j["runs_json"]):
            if (r["run"], stage) in seen:
                continue
            seen.add((r["run"], stage))
            rows.append(("downloaded", r["run"], r.get("bytes", ""), "ENA FASTQ (RNA-seq)", f"cloud/driver.py ({stage})"))

import ast  # noqa: E402
for f in sorted(glob.glob(os.path.join(ROOT, "cloud", "bundle_dna", "jobs", "*.json")) +
                glob.glob(os.path.join(ROOT, "cloud", "bundle_dna", "done_on_mac", "*.json"))):
    d = json.load(open(f))
    rows.append(("downloaded", d["assembly_url"], "", f"assembly for the source-read check {d['job_id']}", "cloud/dna_check.py"))
    for kind in ("illumina", "ont"):
        urls = d.get(kind) or []
        urls = ast.literal_eval(urls) if isinstance(urls, str) else urls
        accs = sorted({os.path.basename(u).split("_")[0].split(".")[0] for u in urls})
        for acc in accs:
            rows.append(("downloaded", acc, d.get(f"{kind}_bytes", "") if len(accs) == 1 else "",
                         f"ENA FASTQ ({kind}), source reads of {d['job_id']}", "cloud/dna_check.py"))
for acc, what in (("GCA_035920445.1", "Nocardia sp. NBC_01730 genome for the misassembly controls"),
                  ("SRR21428551", "Illumina reads for the misassembly controls"),
                  ("SRR21563524", "nanopore reads for the misassembly controls")):
    rows.append(("downloaded", acc, "", what, "controls/dna/build_refs.py, run_controls.sh"))

pf = json.load(open(os.path.join(ROOT, "proteomics", "results", "pride_files_PXD025008.json")))
for x in pf:
    cat = (x["fileCategory"] or {}).get("value")
    if cat in ("RAW", "PEAK", "FASTA", "RESULT"):
        rows.append(("downloaded", "PXD025008/" + x["fileName"], x["fileSizeBytes"], f"PRIDE ({cat})",
                     "proteomics/run_download.sh, run_download_raw.sh"))
for line in open(os.path.join(ROOT, "proteomics", "sycon", "raw", "files.txt")):
    p = line.split()
    if len(p) >= 2 and p[1].isdigit():
        rows.append(("downloaded", "PXD060105/" + p[0], int(p[1]), "PRIDE", "proteomics/sycon/run_download.sh"))

for path, what in (("data/loci/", "locus DNA regions and assemblies fetched by scripts/ (Stage 1, DNA checks)"),
                   ("data/refdb/", "UniProtKB reference set for DIAMOND (Stage 1 orthology), rebuilt by scripts/"),
                   ("data/fasta/*.fasta", "uncompressed copies of the stored .fasta.gz files"),
                   ("cloud_out/", "results synced from the cloud box; a full copy is stored as data/rnaseq_run_results.tar.gz"),
                   ("controls/dna/reads/, controls/dna/work/, controls/dna/refs/*.fa.gz", "control reads, alignments, edited genomes"),
                   ("proteomics/raw/, proteomics/raw_thermo/, proteomics/mzml/", "PXD025008 peak lists, RAW files, converted mzML"),
                   ("proteomics/sycon/raw/*.raw, proteomics/sycon/mzml/", "PXD060105 RAW files and mzML"),
                   ("proteomics/search/out*/, proteomics/sycon/search/out/", "Sage output (results.sage.tsv, 0.1-0.3 GB each)"),
                   ("proteomics/tools/", "Sage and ThermoRawFileParser release binaries (not redistributed)")):
    rows.append(("local-only", path, "", what, "git-ignored"))

with open(os.path.join(ROOT, "DATA_MANIFEST.tsv"), "w") as f:
    f.write("kind\titem\tbytes\tsource\tused_by\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")
st = sum(r[2] for r in rows if r[0] == "stored")
dl = sum(int(r[2]) for r in rows if r[0] == "downloaded" and str(r[2]).isdigit())
print(f"{sum(1 for r in rows if r[0] == 'stored')} stored files, {st / 1e6:.1f} MB; "
      f"{sum(1 for r in rows if r[0] == 'downloaded')} downloaded items, {dl / 1e9:.1f} GB with known size")
