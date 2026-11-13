"""Count and list NCBI Protein entries longer than 35,991 aa via E-utilities.

No API key, no email/tool parameters; <= 2.5 requests/second (see http_util).
NCBI Protein is a superset of what goes into nr (GenPept/INSDC translations,
RefSeq, PDB, Swiss-Prot, PIR, PRF). Category is derived from the accession.

Outputs (data/):
  ncbi_protein_ge35992.tsv   one row per NCBI Protein record (esummary fields)
  ncbi_counts.tsv            esearch counts for several queries
"""
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
from http_util import get  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def esearch(term, retmax=0):
    r = json.loads(get(EU + "esearch.fcgi", {"db": "protein", "term": term,
                                              "retmax": retmax, "retmode": "json"}))
    return int(r["esearchresult"]["count"]), r["esearchresult"].get("idlist", [])


queries = {
    "all_gt_titin_isoform": "35992:100000000[SLEN]",
    "all_gt_pkzilla1": "45213:100000000[SLEN]",
    "refseq_gt_titin_isoform": "35992:100000000[SLEN] AND refseq[filter]",
    "refseq_gt_pkzilla1": "45213:100000000[SLEN] AND refseq[filter]",
    "swissprot_gt_titin_isoform": "35992:100000000[SLEN] AND swissprot[filter]",
    "all_gt_titin_canonical": "34351:100000000[SLEN]",
}
counts = {}
for k, q in queries.items():
    counts[k] = esearch(q)[0]
    print(k, counts[k])
with open(os.path.join(D, "ncbi_counts.tsv"), "w") as f:
    f.write("query_name\tquery\tcount\n")
    for k, q in queries.items():
        f.write(f"{k}\t{q}\t{counts[k]}\n")

n, ids = esearch(queries["all_gt_titin_isoform"], retmax=10000)
assert len(ids) == n, (len(ids), n)
rows = []
for i in range(0, len(ids), 200):
    chunk = ids[i:i + 200]
    r = json.loads(get(EU + "esummary.fcgi", data={"db": "protein", "id": ",".join(chunk),
                                                    "retmode": "json"}))
    for uid in r["result"]["uids"]:
        d = r["result"][uid]
        acc = d.get("accessionversion", "")
        if re.match(r"^(XP|NP|YP|WP|AP|ZP)_", acc):
            cat = "RefSeq-" + acc[:2]
        elif re.match(r"^[A-Z]{3}\d{5,}", acc) or re.match(r"^[A-Z]{1,2}\d{5,}", acc):
            cat = "INSDC"
        elif re.match(r"^\d[A-Z0-9]{3}_", acc):
            cat = "PDB"
        else:
            cat = "other"
        rows.append({
            "uid": uid, "accession": acc, "category": cat, "length": d.get("slen"),
            "taxid": d.get("taxid"), "organism": d.get("organism", ""),
            "title": d.get("title", "").replace("\t", " "),
            "sourcedb": d.get("sourcedb", ""), "createdate": d.get("createdate", ""),
            "updatedate": d.get("updatedate", ""), "status": d.get("status", ""),
            "replacedby": d.get("replacedby", ""), "moltype": d.get("moltype", ""),
            "completeness": d.get("completeness", ""), "genome": d.get("genome", ""),
        })
rows.sort(key=lambda r: -int(r["length"]))
cols = list(rows[0].keys())
with open(os.path.join(D, "ncbi_protein_ge35992.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
print("rows", len(rows))
