"""Positive control outside the > 35,991-aa set: human titin, longest RefSeq isoform
(NCBI Gene 7273, GRCh38). PKZILLA-1 (KAL1508129.1) is already in the main set.
Writes data/loci/control_specs.jsonl, data/fasta/controls.fasta and adds the locus DNA
to data/loci/dna_index.json so 11_locus_dna_checks.py and 12_orthology.py include it.
"""
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from http_util import get  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
LD = os.path.join(D, "loci")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
GENE = "7273"

rep = json.loads(get(f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{GENE}/product_report"))
best = None
for t in rep["reports"][0]["product"]["transcripts"]:
    p = t.get("protein", {})
    if p.get("accession_version", "").startswith("NP_") and t.get("genomic_locations"):
        gl = [g for g in t["genomic_locations"] if g["genomic_accession_version"].startswith("NC_000002")]
        if gl and (best is None or int(p["length"]) > int(best[0]["protein"]["length"])):
            best = (t, gl[0])
t, gl = best
xp = t["protein"]["accession_version"]
ex = sorted(gl["exons"], key=lambda e: int(e["order"]))
cds = t["cds"]["range"][0]
cb, ce = int(cds["begin"]), int(cds["end"])
strand = "+" if gl["genomic_range"]["orientation"] == "plus" else "-"
segs, pos = [], 0
for e in ex:
    b, e_ = int(e["begin"]), int(e["end"])
    L = e_ - b + 1
    lo, hi = max(cb, pos + 1), min(ce, pos + L)
    if lo <= hi:
        if strand == "+":
            segs.append((b + lo - pos - 1, b + hi - pos - 1))
        else:
            segs.append((e_ - (hi - pos - 1), e_ - (lo - pos - 1)))
    pos += L
seq = "".join(get(EU + "efetch.fcgi", {"db": "protein", "id": xp, "rettype": "fasta", "retmode": "text"})
              .decode().split("\n")[1:]).strip()
md5 = "control_human_titin_" + xp
spec = {"md5": md5, "rep_accession": xp, "length": len(seq), "organism": "Homo sapiens",
        "superkingdom": "Eukaryota", "source": f"datasets-gene:{GENE}", "contig": gl["genomic_accession_version"],
        "strand": strand, "segments": segs, "n_contigs_in_model": 1, "note": "positive control (curated RefSeq)"}
open(os.path.join(LD, "control_specs.jsonl"), "w").write(json.dumps(spec) + "\n")
open(os.path.join(D, "fasta", "controls.fasta"), "w").write(f">{md5}\n{seq}\n")
lo, hi = min(a for a, b in segs), max(b for a, b in segs)
a, b = lo - 5000, hi + 5000
lid = f"{spec['contig']}_{a}_{b}"
txt = get(EU + "efetch.fcgi", {"db": "nuccore", "id": spec["contig"], "rettype": "fasta", "retmode": "text",
                               "seq_start": a, "seq_stop": b}, timeout=300).decode()
with gzip.open(os.path.join(LD, "dna", lid + ".fa.gz"), "wt") as f:
    f.write(txt)
idx_path = os.path.join(LD, "dna_index.json")
idx = json.load(open(idx_path))
idx[lid] = {"contig": spec["contig"], "start": a, "end": b, "how": "ncbi", "md5s": [md5]}
json.dump(idx, open(idx_path, "w"))
cl_path = os.path.join(LD, "contig_lengths.json")
cl = json.load(open(cl_path))
cl[spec["contig"]] = 242193529   # GRCh38 chr2 (NC_000002.12) length
json.dump(cl, open(cl_path, "w"))
print(xp, len(seq), "aa", len(segs), "CDS exons", spec["contig"], lo, hi, strand)
