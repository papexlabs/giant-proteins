#!/usr/bin/env python3
"""Per-species genome preparation for stage B (run by driver.py, in the background).

Reads the downloaded assembly (genome.fna.gz), places every locus on it (contig name, or
the 300-bp probe when the assembly names the contig differently, e.g. a CM chromosome
accession for a WGS scaffold), and writes:
  loci_found.json  loci with genome coordinates (relocated if needed)
  loci.bed         locus spans +- 2 kb, used by samtools to keep only locus alignments
Usage: prep_genome.py genome.fna.gz bundle_loci.json out_dir
"""
import gzip
import json
import os
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
genome, loci_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
loci = json.load(open(loci_path))
want = {L["contig"] for L in loci}

names = set()
with gzip.open(genome, "rt") as fi:
    for line in fi:
        if line[0] == ">":
            names.add(line[1:].split()[0])
missing = [L for L in loci if L["contig"] not in names]
if missing:                                   # probe search needs sequences
    seqs, name, buf = {}, None, []
    with gzip.open(genome, "rt") as fi:
        for line in fi:
            if line[0] == ">":
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
        if name:
            seqs[name] = "".join(buf)
    for L in missing:
        probe = L.get("probe", "").upper()
        if not probe:
            continue
        rc = probe.translate(COMP)[::-1]
        for cn, cs in seqs.items():
            i = cs.find(probe)
            if i >= 0:
                off = i + 1 - L["probe_pos"]
                L["segments"] = [[a + off, b + off] for a, b in L["segments"]]
                L["contig_original"], L["contig"] = L["contig"], cn
                break
            j = cs.find(rc)
            if j >= 0:
                off = j + len(probe) + L["probe_pos"] - 1
                L["segments"] = [[off - b, off - a] for a, b in L["segments"]]
                L["strand"] = "-" if L["strand"] == "+" else "+"
                L["contig_original"], L["contig"] = L["contig"], cn
                break
    names = set(seqs)
found = [L for L in loci if L["contig"] in names]
json.dump(found, open(os.path.join(out, "loci_found.json"), "w"))
with open(os.path.join(out, "loci.bed"), "w") as f:
    for L in found:
        lo = min(a for a, b in L["segments"])
        hi = max(b for a, b in L["segments"])
        f.write(f"{L['contig']}\t{max(0, lo - 1 - 2000)}\t{hi + 2000}\t{L['locus_id']}\n")
moved = [f"{L['contig_original']}->{L['contig']}" for L in found if "contig_original" in L]
print(json.dumps({"loci": len(loci), "found": len(found), "relocated": moved}))
