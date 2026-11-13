"""Fetch UniParc and UniProtKB entries longer than the longest human titin isoform.

Outputs (data/):
  raw/uniparc_ge35992.json      raw UniParc JSON stream (sequences, xrefs, features)
  uniparc_ge35992.tsv           one row per UniParc sequence
  uniprotkb_ge35992.tsv         UniProtKB (Swiss-Prot + TrEMBL) entries, UniProt fields
  fasta/uniparc_ge35992.fasta   sequences
"""
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(__file__))
from http_util import get  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
os.makedirs(os.path.join(D, "raw"), exist_ok=True)
os.makedirs(os.path.join(D, "fasta"), exist_ok=True)
MIN_LEN = 35992  # strictly longer than titin's longest isoform (35,991 aa)

# ---- UniParc -------------------------------------------------------------
raw = get("https://rest.uniprot.org/uniparc/stream",
          {"query": f"length:[{MIN_LEN} TO *]", "format": "json"})
open(os.path.join(D, "raw", "uniparc_ge35992.json"), "wb").write(raw)
entries = json.loads(raw)["results"]
# The stream returns "light" entries without cross-references; fetch them per UPI.
XR_CACHE = os.path.join(D, "raw", "uniparc_xrefs.json")
xr_all = json.load(open(XR_CACHE)) if os.path.exists(XR_CACHE) else {}
for i, e in enumerate(entries):
    upi = e["uniParcId"]
    if i % 25 == 0:
        json.dump(xr_all, open(XR_CACHE, "w"))
    if upi not in xr_all:
        xr_all[upi] = json.loads(get(f"https://rest.uniprot.org/uniparc/{upi}/databases",
                                     {"format": "json", "size": 500}))["results"]
json.dump(xr_all, open(XR_CACHE, "w"))

rows = []
with open(os.path.join(D, "fasta", "uniparc_ge35992.fasta"), "w") as fa:
    for e in entries:
        upi = e["uniParcId"]
        seq = e["sequence"]["value"]
        xr = xr_all[upi]
        active = [x for x in xr if x.get("active")]
        orgs = Counter()
        for x in xr:
            o = x.get("organism") or {}
            if o.get("taxonId"):
                orgs[(o.get("taxonId"), o.get("scientificName", ""))] += 1
        (taxid, org), _ = orgs.most_common(1)[0] if orgs else ((None, ""), 0)
        names = [x.get("proteinName") for x in xr if x.get("proteinName")]
        genes = [x.get("geneName") for x in xr if x.get("geneName")]
        pfam = Counter(f["databaseId"] for f in e.get("sequenceFeatures", []) if f["database"] == "Pfam")
        pfam_names = {}
        for f in e.get("sequenceFeatures", []):
            if f["database"] == "Pfam" and f.get("interproGroup"):
                pfam_names[f["databaseId"]] = f["interproGroup"]["name"]
        rows.append({
            "upi": upi,
            "length": e["sequence"]["length"],
            "organism": org,
            "taxid": taxid or "",
            "n_organisms": len(orgs),
            "active": "yes" if active else "no",
            "active_dbs": ";".join(sorted({x["database"] for x in active})),
            "active_ids": ";".join(f'{x["database"]}:{x["id"]}' for x in active),
            "inactive_dbs": ";".join(sorted({x["database"] for x in xr if not x.get("active")})),
            "inactive_ids": ";".join(f'{x["database"]}:{x["id"]}' for x in xr if not x.get("active"))[:500],
            "protein_name": (Counter(names).most_common(1)[0][0] if names else ""),
            "gene": (Counter(genes).most_common(1)[0][0] if genes else ""),
            "first_seen": e.get("oldestCrossRefCreated", ""),
            "last_updated": e.get("mostRecentCrossRefUpdated", ""),
            "pfam_top": ";".join(f"{k}({pfam_names.get(k,'?')})x{v}" for k, v in pfam.most_common(6)),
        })
        fa.write(f">{upi} len={len(seq)} org={org.replace(' ', '_')}\n")
        for i in range(0, len(seq), 80):
            fa.write(seq[i:i + 80] + "\n")

rows.sort(key=lambda r: -r["length"])
cols = list(rows[0].keys())
with open(os.path.join(D, "uniparc_ge35992.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
print("uniparc entries:", len(rows))

# ---- UniProtKB ---------------------------------------------------------------
fields = ",".join([
    "accession", "id", "reviewed", "length", "organism_name", "organism_id", "lineage",
    "protein_name", "gene_names", "protein_existence", "annotation_score", "fragment",
    "xref_embl", "xref_refseq", "xref_proteomes", "xref_pfam", "xref_interpro",
    "date_created", "date_sequence_modified", "version", "uniparc_id",
])
tsv = get("https://rest.uniprot.org/uniprotkb/stream",
          {"query": f"length:[{MIN_LEN} TO *]", "format": "tsv", "fields": fields})
open(os.path.join(D, "uniprotkb_ge35992.tsv"), "wb").write(tsv)
print("uniprotkb rows:", tsv.decode().count("\n") - 1)
