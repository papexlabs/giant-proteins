"""Merge NCBI Protein and UniParc giant sequences into one non-redundant table
and compute cheap automatic flags per unique sequence.

Unique key = MD5 of the uppercase sequence without terminal '*'.
Flags (all computed here, no external tools):
  x_frac            fraction of X residues (gaps/ambiguous codons in the source)
  internal_stop     '*' inside the sequence
  starts_M          first residue is Met
  lowcomp_frac      fraction of residues inside 12-aa windows with Shannon entropy < 2.2 bits
                    (SEG-like trigger threshold; an approximation, not SEG itself)
  rep25_frac        fraction of residues covered by 25-mers that occur >= 2 times in the
                    same sequence (near-identical internal duplication: real tandem repeats
                    or repeat expansion / duplicated exon chains in the model)
  pfam_*            from UniParc sequence features (InterPro member database matches
                    pre-computed by EBI), when the sequence is in UniParc
  arch_class        titin-like / PKS / NRPS / PKS+NRPS / other / none (rule based)
Outputs:
  data/giant_unique_sequences.tsv
"""
import csv
import hashlib
import json
import math
import os
import re
from collections import Counter, defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")


def read_fasta(p):
    name, buf = None, []
    for line in open(p):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def md5(s):
    return hashlib.md5(s.upper().rstrip("*").encode()).hexdigest()


def lowcomp_frac(s, w=12, thr=2.2):
    n = len(s)
    if n < w:
        return 0.0
    mask = bytearray(n)
    c = Counter(s[:w])

    def ent(c):
        return -sum(v / w * math.log2(v / w) for v in c.values() if v)
    for i in range(n - w + 1):
        if i:
            c[s[i - 1]] -= 1
            c[s[i + w - 1]] += 1
        if ent(c) < thr:
            mask[i:i + w] = b"\x01" * w
    return sum(mask) / n


def rep_frac(s, k=25):
    seen = Counter(s[i:i + k] for i in range(len(s) - k + 1))
    mask = bytearray(len(s))
    for i in range(len(s) - k + 1):
        if seen[s[i:i + k]] > 1:
            mask[i:i + k] = b"\x01" * k
    return sum(mask) / max(1, len(s))


TITIN = {"PF07679", "PF00041", "PF13927", "PF00047", "PF13895", "PF00069", "PF07714"}
PKS = {"PF00109", "PF02801", "PF00698", "PF08659", "PF14765", "PF08240", "PF13602"}
NRPS = {"PF00501", "PF00668", "PF13193"}

# --- sources --------------------------------------------------------------------
ncbi_meta = {r["accession"]: r for r in csv.DictReader(open(os.path.join(D, "ncbi_records_parsed.tsv")), delimiter="\t")}
ncbi_sum = {r["accession"]: r for r in csv.DictReader(open(os.path.join(D, "ncbi_protein_ge35992.tsv")), delimiter="\t")}
up_rows = {r["upi"]: r for r in csv.DictReader(open(os.path.join(D, "uniparc_ge35992.tsv")), delimiter="\t")}
ukb = {r["Entry"]: r for r in csv.DictReader(open(os.path.join(D, "uniprotkb_ge35992.tsv")), delimiter="\t")}
up_raw = {e["uniParcId"]: e for e in json.load(open(os.path.join(D, "raw", "uniparc_ge35992.json")))["results"]}

seqs = {}
members = defaultdict(lambda: {"ncbi": [], "upi": []})
for acc, s in read_fasta(os.path.join(D, "fasta", "ncbi_ge35992.fasta")):
    h = md5(s)
    seqs[h] = s.rstrip("*")
    members[h]["ncbi"].append(acc)
for upi, s in read_fasta(os.path.join(D, "fasta", "uniparc_ge35992.fasta")):
    h = md5(s)
    seqs.setdefault(h, s)
    members[h]["upi"].append(upi)

ukb_by_upi = defaultdict(list)
for acc, r in ukb.items():
    ukb_by_upi[r["UniParc"].strip(";")].append(acc)

out = []
for h, s in seqs.items():
    m = members[h]
    upi = m["upi"][0] if m["upi"] else ""
    ur = up_rows.get(upi, {})
    ncbi = sorted(m["ncbi"], key=lambda a: (not re.match(r"^(XP|NP|WP)_", a), a))
    rep = ncbi[0] if ncbi else (ukb_by_upi.get(upi, [""])[0] or upi)
    nm = ncbi_meta.get(ncbi[0], {}) if ncbi else {}
    ns = ncbi_sum.get(ncbi[0], {}) if ncbi else {}
    ukb_accs = ukb_by_upi.get(upi, [])
    uk = ukb.get(ukb_accs[0], {}) if ukb_accs else {}
    # Pfam coverage from UniParc features
    feats = up_raw.get(upi, {}).get("sequenceFeatures", []) if upi else []
    pf = [f for f in feats if f["database"] == "Pfam"]
    cov = bytearray(len(s))
    fams = Counter()
    for f in pf:
        fams[f["databaseId"]] += len(f["locations"])
        for loc in f["locations"]:
            cov[loc["start"] - 1:loc["end"]] = b"\x01" * (loc["end"] - loc["start"] + 1)
    # longest stretch without any InterPro member-db feature (Pfam/SMART/Gene3D/...)
    anycov = bytearray(len(s))
    for f in feats:
        if f["database"] in ("MobiDBLite", "Coils"):
            continue
        for loc in f["locations"]:
            anycov[loc["start"] - 1:loc["end"]] = b"\x01" * (loc["end"] - loc["start"] + 1)
    longest_gap = max((len(x) for x in re.split(b"\x01+", bytes(anycov))), default=0) if feats else ""
    # ordered Pfam architecture, consecutive repeats collapsed (e.g. I-set*12)
    hits = sorted((loc["start"], f["databaseId"]) for f in pf for loc in f["locations"])
    arch_runs = []
    for _, fam in hits:
        if arch_runs and arch_runs[-1][0] == fam:
            arch_runs[-1][1] += 1
        else:
            arch_runs.append([fam, 1])
    arch_string = ">".join(f"{a}*{n}" if n > 1 else a for a, n in arch_runs)
    fs = set(fams)
    if not upi or not feats:
        arch = "not_in_uniparc" if not upi else "no_features"
    elif fs & PKS and fs & NRPS:
        arch = "PKS+NRPS"
    elif len(fs & PKS) >= 2:
        arch = "PKS"
    elif fs & NRPS:
        arch = "NRPS"
    elif fs & {"PF07679", "PF00041", "PF13927"} and sum(fams[x] for x in fs & TITIN) >= 50:
        arch = "titin-like"
    elif fams:
        arch = "other"
    else:
        arch = "no_pfam"
    title = ns.get("title", "") or ur.get("protein_name", "")
    org = nm.get("organism") or ns.get("organism") or ur.get("organism", "")
    patent = (nm.get("patent") == "yes") or "patent" in title.lower() or ur.get("organism", "") == "synthetic construct"
    live_ncbi = bool(ncbi)
    active_up = ur.get("active") == "yes"
    out.append({
        "md5": h, "rep_accession": rep, "length": len(s), "organism": org,
        "superkingdom": nm.get("superkingdom", ""), "title": title,
        "ncbi_accessions": ";".join(ncbi), "upi": ";".join(m["upi"]),
        "uniprotkb": ";".join(ukb_accs),
        "uniprotkb_PE": uk.get("Protein existence", ""), "uniprotkb_reviewed": uk.get("Reviewed", ""),
        "in_ncbi_protein": "yes" if live_ncbi else "no",
        "uniparc_active": ur.get("active", ""), "uniparc_active_dbs": ur.get("active_dbs", ""),
        "uniparc_inactive_dbs": ur.get("inactive_dbs", ""),
        "currently_live": "yes" if (live_ncbi or active_up) else "no",
        "patent_or_synthetic": "yes" if patent else "no",
        "refseq_low_quality": nm.get("low_quality", ""),
        "partial_in_title": "yes" if "partial" in title.lower() else "no",
        "MAG": "yes" if "MAG:" in title else "no",
        "n_cds_parts": nm.get("n_cds_parts", ""), "genomic_span_bp": nm.get("genomic_span_bp", ""),
        "locus_contigs": nm.get("locus_contigs", ""), "gene_exoncount": nm.get("gene_exoncount", ""),
        "rnaseq_feature_coverage_pct": nm.get("rnaseq_feature_coverage_pct", ""),
        "rnaseq_samples_all_introns": nm.get("rnaseq_samples_all_introns", ""),
        "refseq_model_evidence": nm.get("refseq_model_evidence", ""),
        "x_frac": round(s.count("X") / len(s), 5), "internal_stop": "yes" if "*" in s else "no",
        "starts_M": "yes" if s[:1] == "M" else "no",
        "lowcomp_frac": round(lowcomp_frac(s), 4), "rep25_frac": round(rep_frac(s), 4),
        "pfam_cov_frac": round(sum(cov) / len(s), 4) if feats else "",
        "longest_unannotated_aa": longest_gap,
        "n_pfam_families": len(fams) if feats else "",
        "pfam_top": ";".join(f"{k}x{v}" for k, v in fams.most_common(8)),
        "arch_class": arch,
        "pfam_arch_ordered": arch_string[:400],
    })

# lineage for sequences that have no NCBI record (UniParc-only): UniProt taxonomy API
import sys  # noqa: E402
sys.path.insert(0, os.path.dirname(__file__))
from http_util import get  # noqa: E402
TAXC = os.path.join(D, "raw", "taxonomy_cache.json")
taxc = json.load(open(TAXC)) if os.path.exists(TAXC) else {}
for r in out:
    if r["superkingdom"] or not r["upi"]:
        continue
    tid = str(up_rows[r["upi"].split(";")[0]]["taxid"])
    if not tid:
        continue
    if tid not in taxc:
        try:
            t = json.loads(get(f"https://rest.uniprot.org/taxonomy/{tid}", {"format": "json"}))
            taxc[tid] = [x["scientificName"] for x in t.get("lineage", [])]
        except Exception:  # noqa: BLE001
            taxc[tid] = []
    lin = taxc[tid]
    r["superkingdom"] = next((x for x in lin if x in ("Eukaryota", "Bacteria", "Archaea", "Viruses")),
                             "other/unknown")
json.dump(taxc, open(TAXC, "w"))

out.sort(key=lambda r: -r["length"])
with open(os.path.join(D, "fasta", "unique_ge35992.fasta"), "w") as fa:
    for r in out:
        s = seqs[r["md5"]]
        fa.write(f">{r['md5']} {r['rep_accession']} len={r['length']}\n")
        fa.write("\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
cols = list(out[0].keys())
with open(os.path.join(D, "giant_unique_sequences.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in out:
        f.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
print("unique sequences:", len(out), "live:", sum(r["currently_live"] == "yes" for r in out))
