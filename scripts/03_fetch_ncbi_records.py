"""Fetch GenPept flat files for all NCBI Protein records > 35,991 aa and parse
locus-level metadata (coded_by -> exon count and genomic span, lineage, RefSeq
model evidence for XP_ records via the XM_ feature table and NCBI Gene).

Outputs:
  data/raw/ncbi_gp.txt.gz          concatenated GenPept records
  data/raw/xm_ft.txt.gz            feature tables of the XM_ mRNAs (no sequence)
  data/fasta/ncbi_ge35992.fasta    sequences
  data/ncbi_records_parsed.tsv     one row per record
"""
import csv
import gzip
import hashlib
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
from http_util import get  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

rows = list(csv.DictReader(open(os.path.join(D, "ncbi_protein_ge35992.tsv")), delimiter="\t"))
uids = [r["uid"] for r in rows]

GP = os.path.join(D, "raw", "ncbi_gp.txt.gz")
if not os.path.exists(GP):
    part_dir = os.path.join(D, "raw", "gp_parts")
    os.makedirs(part_dir, exist_ok=True)
    for i in range(0, len(uids), 5):
        pp = os.path.join(part_dir, f"{i:05d}.gp")
        if os.path.exists(pp):
            continue
        t = get(EU + "efetch.fcgi", data={"db": "protein", "id": ",".join(uids[i:i + 5]),
                                          "rettype": "gp", "retmode": "text"}, timeout=60).decode()
        open(pp, "w").write(t)
        print("gp batch", i, len(t), flush=True)
    with gzip.open(GP, "wt") as f:
        for fn in sorted(os.listdir(part_dir)):
            f.write(open(os.path.join(part_dir, fn)).read())
text = gzip.open(GP, "rt").read()
records = [r for r in text.split("\n//\n") if r.strip().startswith("LOCUS")]


def section(rec, name):
    """Return the text block of a top-level section (continuation lines included)."""
    m = re.search(rf"^{name}\s+(.*?)(?=^\S)", rec, re.S | re.M)
    return re.sub(r"\s+", " ", m.group(1)).strip() if m else ""


def qualifiers(rec, key):
    return [re.sub(r"\s+", "", m) if key == "coded_by" else re.sub(r"\s+", " ", m)
            for m in re.findall(rf'/{key}="(.*?)"', rec, re.S)]


def span_of(loc):
    parts = re.findall(r"([A-Z]{1,6}_?\d+(?:\.\d+)?):<?(\d+)\.\.>?(\d+)", loc)
    if not parts:
        return "", 0, 0, 0
    contigs = sorted({p[0] for p in parts})
    coords = [int(x) for p in parts for x in p[1:]]
    return ";".join(contigs), len(parts), min(coords), max(coords)


out = []
fa = open(os.path.join(D, "fasta", "ncbi_ge35992.fasta"), "w")
for rec in records:
    acc = section(rec, "VERSION").split()[0]
    org_block = re.search(r"^  ORGANISM\s+(.*?)(?=^REFERENCE|^COMMENT|^FEATURES)", rec, re.S | re.M)
    org_lines = org_block.group(1).split("\n") if org_block else [""]
    organism = org_lines[0].strip()
    lineage = re.sub(r"\s+", " ", " ".join(org_lines[1:])).strip().rstrip(".")
    seq = re.sub(r"[^A-Za-z*]", "", rec.split("\nORIGIN", 1)[1]) .upper() if "\nORIGIN" in rec else ""
    coded = qualifiers(rec, "coded_by")
    contigs, nparts, start, end = span_of(coded[0]) if coded else ("", 0, 0, 0)
    xrefs = qualifiers(rec, "db_xref")
    gene_id = next((x.split(":")[1] for x in xrefs if x.startswith("GeneID:")), "")
    comment = section(rec, "COMMENT")
    notes = " | ".join(n for n in qualifiers(rec, "note") if "[structural motif]" not in n
                       and not re.search(r"; c[dl]\d+$", n))[:600]
    flat = re.sub(r"\s+", " ", rec.split("\nORIGIN", 1)[0])
    mm = re.search(r"modified relative to its source genomic sequence to represent the inferred CDS: ([^\"]*)", flat)
    out.append({
        "accession": acc, "length": len(seq), "organism": organism,
        "superkingdom": lineage.split(";")[0].strip() if lineage else "",
        "lineage": lineage, "definition": section(rec, "DEFINITION"),
        "dbsource": section(rec, "DBSOURCE")[:200], "keywords": section(rec, "KEYWORDS"),
        "coded_by": (coded[0][:300] if coded else ""), "n_cds_parts": nparts,
        "locus_contigs": contigs, "locus_start": start, "locus_end": end,
        "genomic_span_bp": (end - start + 1) if nparts else "",
        "gene_id": gene_id,
        "low_quality": "yes" if "LOW QUALITY" in rec[:3000] else "no",
        "patent": "yes" if "patent" in section(rec, "DEFINITION").lower() or "PAT" in section(rec, "LOCUS") else "no",
        "comment_head": comment[:300], "notes": notes,
        "model_modifications": mm.group(1).strip() if mm else "",
        "md5": hashlib.md5(seq.replace("*", "").encode()).hexdigest(),
    })
    fa.write(f">{acc} len={len(seq)}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
fa.close()

# RefSeq model evidence for XP_ records: the XM_ GenBank record (fetched truncated
# to 1 bp so no sequence is downloaded) keeps the COMMENT and the mRNA /note with
# "Supporting evidence includes similarity to: ... RNAseq alignments ...".
XM = os.path.join(D, "raw", "xm_gb_trunc.txt.gz")
xms = sorted({o["coded_by"].split(":")[0] for o in out if o["coded_by"].startswith("XM_")})
if not os.path.exists(XM) and xms:
    chunks = []
    for i in range(0, len(xms), 20):
        chunks.append(get(EU + "efetch.fcgi", data={"db": "nuccore", "id": ",".join(xms[i:i + 20]),
                                                     "rettype": "gb", "retmode": "text",
                                                     "seq_start": 1, "seq_stop": 1}, timeout=60).decode())
    with gzip.open(XM, "wt") as f:
        f.write("".join(chunks))
evidence, xm_attr = {}, {}
if xms:
    for rec in gzip.open(XM, "rt").read().split("\n//\n"):
        m = re.search(r"^VERSION\s+(XM_\d+\.\d+)", rec, re.M)
        if not m:
            continue
        flat = re.sub(r"\s+", " ", rec)
        ev = re.search(r"Supporting evidence includes similarity to: ([^\"]*)", flat)
        evidence[m.group(1)] = ev.group(1).strip() if ev else ""
        attrs = re.search(r"##RefSeq-Attributes-START## (.*?) ##RefSeq-Attributes-END##", flat)
        xm_attr[m.group(1)] = attrs.group(1).strip() if attrs else ""

# NCBI Gene: genomic location and exon count for records with a GeneID.
gids = sorted({o["gene_id"] for o in out if o["gene_id"]})
gene = {}
for i in range(0, len(gids), 200):
    r = json.loads(get(EU + "esummary.fcgi", data={"db": "gene", "id": ",".join(gids[i:i + 200]),
                                                    "retmode": "json"}))
    for uid in r["result"].get("uids", []):
        g = r["result"][uid]
        gi = (g.get("genomicinfo") or [{}])[0]
        gene[uid] = {"gene_symbol": g.get("name", ""), "gene_desc": g.get("description", ""),
                     "gene_chr": gi.get("chraccver", ""), "gene_start": gi.get("chrstart", ""),
                     "gene_stop": gi.get("chrstop", ""), "gene_exoncount": gi.get("exoncount", "")}

for o in out:
    xm = o["coded_by"].split(":")[0] if o["coded_by"].startswith("XM_") else ""
    o["refseq_model_evidence"] = evidence.get(xm, "")
    o["refseq_attributes"] = xm_attr.get(xm, "")
    m = re.search(r"substituted (\d+) bases? at (\d+) genomic stop codons?", o["model_modifications"])
    o["genomic_stops_substituted"] = m.group(2) if m else "0"
    m = re.search(r"(\d+)% coverage of the annotated genomic feature by RNAseq alignments", o["refseq_model_evidence"])
    o["rnaseq_feature_coverage_pct"] = m.group(1) if m else ""
    m = re.search(r"including (\d+) samples? with support for all annotated introns", o["refseq_model_evidence"])
    o["rnaseq_samples_all_introns"] = m.group(1) if m else ("0" if o["rnaseq_feature_coverage_pct"] else "")
    o.update(gene.get(o["gene_id"], {"gene_symbol": "", "gene_desc": "", "gene_chr": "", "gene_start": "",
                                      "gene_stop": "", "gene_exoncount": ""}))

cols = list(out[0].keys())
with open(os.path.join(D, "ncbi_records_parsed.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for o in sorted(out, key=lambda o: -o["length"]):
        f.write("\t".join(str(o[c]).replace("\t", " ") for c in cols) + "\n")
print("records", len(out), "xm evidence", sum(bool(v) for v in evidence.values()), "genes", len(gene))
