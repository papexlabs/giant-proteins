"""Stage 0 table: one row per unique giant sequence (> 35,991 aa) with database
membership, metadata and cheap automatic flags, plus a provisional triage label.

The triage label is NOT an evidence grade. It only sorts entries for Stage 1:
  T-excl      synthetic/patent, no live record, or marked partial
  T-repeat    sequence pattern of a translated tandem DNA repeat or low-complexity
              run with (almost) no recognised domains
  T-arch      recognised multi-domain architecture covering >= 30 % of the
              length (PKS/NRPS, Ig/FN3 titin-like, adhesin repeats, dynein)
              and no repeat-translation pattern
  T-unclear   everything else
Outputs: data/stage0_table.tsv, data/stage0_loci.tsv
"""
import csv
import os
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")


def load(name, key):
    return {r[key]: r for r in csv.DictReader(open(os.path.join(D, name)), delimiter="\t")}


uniq = list(csv.DictReader(open(os.path.join(D, "giant_unique_sequences.tsv")), delimiter="\t"))
dom = load("domain_panel_summary.tsv", "md5")
per = load("period_flags.tsv", "md5")
rec = load("ncbi_records_parsed.tsv", "accession")

rows = []
for u in uniq:
    d, p = dom.get(u["md5"], {}), per.get(u["md5"], {})
    accs = [a for a in u["ncbi_accessions"].split(";") if a]
    r0 = rec.get(accs[0], {}) if accs else {}
    panel_cov = float(d.get("panel_cov_frac") or 0)
    pf_cov = float(u["pfam_cov_frac"] or 0) if u["pfam_cov_frac"] else None
    dom_cov = max(panel_cov, pf_cov or 0)
    per_id = float(p.get("top_period_identity") or 0)
    rep25 = float(u["rep25_frac"])
    low = float(u["lowcomp_frac"])
    flags = []
    if u["patent_or_synthetic"] == "yes":
        flags.append("synthetic/patent")
    if u["currently_live"] != "yes":
        flags.append("no-live-record")
    if u["partial_in_title"] == "yes":
        flags.append("partial")
    if (per_id >= 0.6 or rep25 >= 0.8) and dom_cov < 0.05:
        flags.append(f"tandem-repeat-translation(period={p.get('top_period_aa')}aa,id={per_id:.2f})")
    if low >= 0.3:
        flags.append(f"low-complexity({low:.2f})")
    if float(u["x_frac"]) > 0:
        flags.append(f"X-residues({u['x_frac']})")
    if u["internal_stop"] == "yes":
        flags.append("internal-stop")
    if u["starts_M"] != "yes":
        flags.append("no-initial-Met")
    if u["refseq_low_quality"] == "yes":
        flags.append("RefSeq-LOW-QUALITY")
    if r0.get("genomic_stops_substituted", "0") not in ("0", ""):
        flags.append(f"genomic-stop-substituted({r0['genomic_stops_substituted']})")
    if "gap" in r0.get("refseq_attributes", "") or "not found in genome assembly" in r0.get("model_modifications", ""):
        flags.append("model-patches-assembly-gap")
    if "ab initio" in r0.get("refseq_attributes", ""):
        flags.append("RefSeq-ab-initio-part:" + r0["refseq_attributes"].split("::")[-1].strip()[:20])
    if u["MAG"] == "yes":
        flags.append("MAG")
    if dom_cov == 0 and not d.get("panel_hits") and not u["pfam_top"]:
        flags.append("no-domain-detected")

    if {"synthetic/patent", "no-live-record", "partial"} & set(flags):
        tri = "T-excl"
    elif any(f.startswith("tandem-repeat-translation") for f in flags) or (low >= 0.5 and dom_cov < 0.05):
        tri = "T-repeat"
    elif d.get("panel_class", "").split(" ")[0] in ("PKS", "NRPS", "PKS+NRPS", "Ig/FN3-repeat", "dynein-HC-like",
                                                   "bacterial-adhesin-repeat") and dom_cov >= 0.30:
        tri = "T-arch"
    else:
        tri = "T-unclear"

    # locus key: same gene (RefSeq GeneID) or same contig region; otherwise the sequence itself
    if r0.get("gene_id"):
        locus = f"GeneID:{r0['gene_id']}"
    elif r0.get("locus_contigs") and r0.get("locus_start"):
        locus = f"{r0['locus_contigs']}:{r0['locus_start']}-{r0['locus_end']}"
    else:
        locus = f"seq:{u['upi'] or u['md5']}"
    rows.append({
        "rep_accession": u["rep_accession"], "length": u["length"], "organism": u["organism"],
        "superkingdom": u["superkingdom"], "title": u["title"][:90],
        "databases": ";".join(x for x in [
            "NCBI-Protein" if u["in_ncbi_protein"] == "yes" else "",
            "RefSeq" if any(a[:3] in ("XP_", "NP_", "WP_") for a in accs) else "",
            "UniProtKB/TrEMBL" if u["uniprotkb"] else "",
            "UniParc(active)" if u["uniparc_active"] == "yes" else ("UniParc(inactive)" if u["upi"] else ""),
        ] if x),
        "n_ncbi_accessions": len(accs), "ncbi_accessions": u["ncbi_accessions"][:200],
        "uniprotkb": u["uniprotkb"], "uniprotkb_PE": u["uniprotkb_PE"], "upi": u["upi"],
        "live": u["currently_live"], "locus_key": locus,
        "n_cds_parts": u["n_cds_parts"], "gene_exoncount": u["gene_exoncount"],
        "genomic_span_bp": u["genomic_span_bp"],
        "refseq_model_evidence": u["refseq_model_evidence"][:160],
        "rnaseq_feature_coverage_pct": u["rnaseq_feature_coverage_pct"],
        "rnaseq_samples_all_introns": u["rnaseq_samples_all_introns"],
        "x_frac": u["x_frac"], "lowcomp_frac": u["lowcomp_frac"], "rep25_frac": u["rep25_frac"],
        "top_period_aa": p.get("top_period_aa", ""), "top_period_identity": p.get("top_period_identity", ""),
        "panel_class": d.get("panel_class", ""), "panel_cov_frac": d.get("panel_cov_frac", ""),
        "uniparc_pfam_cov_frac": u["pfam_cov_frac"], "uniparc_pfam_top": u["pfam_top"][:120],
        "flags": ";".join(flags), "triage": tri, "md5": u["md5"],
    })

rows.sort(key=lambda r: -int(r["length"]))
cols = list(rows[0].keys())
with open(os.path.join(D, "stage0_table.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")

# loci: collapse isoforms / identical regions; keep the longest sequence per locus
loci = defaultdict(list)
for r in rows:
    loci[r["locus_key"]].append(r)
with open(os.path.join(D, "stage0_loci.tsv"), "w") as f:
    f.write("locus_key\tn_sequences\tmax_length\trep_accession\torganism\ttriage\tpanel_class\tflags\n")
    for k, rs in sorted(loci.items(), key=lambda kv: -int(kv[1][0]["length"])):
        top = rs[0]
        f.write(f"{k}\t{len(rs)}\t{top['length']}\t{top['rep_accession']}\t{top['organism']}\t"
                f"{top['triage']}\t{top['panel_class']}\t{top['flags']}\n")
print("rows", len(rows), "loci", len(loci))
