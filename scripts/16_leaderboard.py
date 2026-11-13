"""Leaderboard of the largest protein entries with their evidence grade (Stage 2).

One row per locus (longest live model per gene or contig region), largest first, from
data/locus_checks.tsv (graded loci) plus the Stage 0 exclusions (grade X: synthetic or
patent, partial, no live record). Columns: rank (A-D only), accession, organism, length,
grade, best grade among the locus' other models, one-line evidence, what is still missing.
The same gene recorded twice in one genome (haplotypes, GenBank + RefSeq, isoforms) is merged
into one row (same_gene_other_records; >= 90 % shared 10-mers, A-C rows only).
PRIDE availability per species comes from data/species_data_availability.tsv (keyword
match on the PRIDE Archive). One dataset was reanalysed (PXD025008, for WP_128698900.1 and its
GenBank twin); its result comes from data/followups/peptide_evidence.tsv.
Output: data/leaderboard.tsv
"""
import csv
import gzip
import os
import re

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
checks = list(csv.DictReader(open(os.path.join(D, "locus_checks.tsv")), delimiter="\t"))
s0 = {r["md5"]: r for r in csv.DictReader(open(os.path.join(D, "stage0_table.tsv")), delimiter="\t")}
rna = {r["md5"]: r for r in csv.DictReader(open(os.path.join(D, "rnaseq_support.tsv")), delimiter="\t")}
avail = {r["organism"]: r for r in csv.DictReader(open(os.path.join(D, "species_data_availability.tsv")),
                                                  delimiter="\t")}
_pp = os.path.join(D, "followups", "peptide_evidence.tsv")
PEPTIDES = {r["accession"]: r for r in csv.DictReader(open(_pp), delimiter="\t")} if os.path.exists(_pp) else {}
ORDER = "ABCDX"


def evidence(r):
    g, reason = r["grade"], r["grade_reason_stage1"]
    v = rna.get(r["md5"], {})
    parts = []
    if g == "D":
        m = re.search(r"(D[1-5] [^;]*)", reason)
        parts.append(m.group(1)[:110] if m else "contradiction")
        if "fusion boundary" in r["grade_change_reason"]:
            fb = r["grade_change_reason"].split("fusion boundary: ")[1]
            parts.append("RNA-seq: " + fb[:120])
        return "; ".join(parts)
    if r["n_exons"] == "1":
        parts.append(f"single ORF {int(r['cds_bp'])//1000} kb, no N, coding signal")
    elif r["n_exons"]:
        parts.append(f"{r['n_exons']} exons, no N in CDS")
    for tag in ("S1 ", "S2 ", "S3 ", "S4 "):
        m = re.search(re.escape(tag) + r"([^;]*)", reason)
        if m:
            parts.append(tag.strip() + " " + m.group(1)[:70])
    dna = cluster_dna(r)
    if len(dna) == 1:
        parts.append("source reads: " + dna[0][1][:110])
    elif dna:
        parts.append("source reads per assembly: " + "; ".join(f"{j}: {v[:70]}" for j, v in dna))
    if v.get("runs"):
        parts.append(f"RNA-seq {v['introns_counted_supported']}/{v['n_introns']} introns "
                     f"({v['introns_unique_support']} unique), {v['verdict']}")
    pe = PEPTIDES.get(r["rep_accession"])
    if pe:
        parts.append(f"peptides: {pe['unique_peptides']} unique at 1 % FDR, all in aa {pe['first_peptide_start']}-"
                     f"{pe['last_peptide_end']} ({pe['dataset']} re-search)")
    if g == "C":
        m = re.search(r"(gene ends < 1 kb[^;]*|weak coding signal[^;]*|locus not resolved[^;]*)", reason)
        if m:
            parts.append(m.group(1)[:80])
    return "; ".join(parts)[:420]


def missing(r):
    g = r["grade"]
    v = rna.get(r["md5"], {})
    euk = r["superkingdom"] == "Eukaryota"
    out = []
    if g == "D":
        return "a corrected gene model or assembly; the current model is contradicted"
    if euk:
        verdict = v.get("verdict", "")
        if verdict not in ("full",):
            if verdict.startswith("no-rnaseq"):
                out.append("RNA-seq (none usable in ENA)")
            elif verdict.startswith("not-tested"):
                out.append("RNA-seq test (" + verdict.split("(")[-1].rstrip(")") + ")")
            elif verdict:
                out.append(f"full-length transcript support (RNA-seq: {verdict})")
            else:
                out.append("RNA-seq test")
        out.append("long-read transcript spanning the whole gene")
    else:
        vs = [v for _, v in cluster_dna(r)]
        if any(v.startswith("supported by the reads") for v in vs):
            pass
        elif any(v.startswith("supported outside") for v in vs):
            kb = min(float(re.search(r"outside ([0-9.]+) kb", v).group(1)) for v in vs if v.startswith("supported outside"))
            out.append(f"long reads or long-insert pairs across the {kb} kb near-identical repeat that the "
                       "short pairs cannot cross (no public long reads for these samples)")
        elif vs:
            out.append("read-level validation of the assembly (" + "; ".join(v.split(":")[0] for v in vs) + ")")
        else:
            out.append("read-level validation of the assembly (source genome/metagenome reads)")
        if not PEPTIDES.get(r["rep_accession"]):
            out.append("any expression data (none public for this strain/MAG)")
    if "S1 " not in r["grade_reason_stage1"]:
        out.append("orthologue of consistent length")
    if "S2 " not in r["grade_reason_stage1"]:
        out.append("independent assembly with the same model")
    ti = r.get("translation_identity", "")
    if ti and float(ti) < 0.95:
        out.append(f"check that the assembly encodes the protein (exon chain translation identity {ti}: "
                   f"codon_start, RefSeq indel correction or gap patch)")
    pe = PEPTIDES.get(r["rep_accession"])
    if pe:
        out.append(f"peptides after aa {pe['last_peptide_end']} of {r['length']} ({pe['dataset']} re-search found none)")
        return "; ".join(out)
    pr = avail.get(r["organism"], {}).get("pride_projects_name_match", "")
    out.append("peptides" + (f" (PRIDE: {pr} project(s) mention the species, not reanalysed)"
                             if pr not in ("", "0", "error") else " (no PRIDE project found for the species)"))
    return "; ".join(out)


def cluster_dna(r):
    """DNA-check verdicts of this model and its near-identical copies in other assemblies."""
    out = []
    for a in copies_of.get(find(r["md5"]), [r["rep_accession"]]):
        x = by_acc.get(a, {})
        if x.get("dna_check_verdict"):
            j = x.get("dna_check_job", "")
            out.append((j.split("_")[1] if j.count("_") >= 2 else j, x["dna_check_verdict"]))
    return sorted(out)


# near-identical copies of one gene in several assemblies of the same species (S2 copies)
# are one leaderboard row: union-find over the copy lists
acc2md5 = {r["rep_accession"]: r["md5"] for r in checks}
by_acc = {r["rep_accession"]: r for r in checks}
parent = {}


def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


for r in checks:
    m = re.search(r"near-identical copy in another assembly: ([^;]*(?:;[A-Za-z0-9_.]+)*)", r["grade_reason_stage1"])
    if m:
        for c in m.group(1).split(";"):
            if c.strip() in acc2md5:
                parent[find(acc2md5[c.strip()])] = find(r["md5"])
copies_of = {}
for r in checks:
    copies_of.setdefault(find(r["md5"]), []).append(r["rep_accession"])

# loci: longest model per locus_key (or copy cluster), plus the best grade among its models
best = {}
loci = {}
for r in checks:
    if r["md5"].startswith("control_"):
        continue
    k = s0[r["md5"]]["locus_key"]
    if len(copies_of[find(r["md5"])]) > 1:
        k = "copies:" + find(r["md5"])
    best[k] = min(best.get(k, "X"), r["grade"], key=ORDER.index)
    if k not in loci or int(r["length"]) > int(loci[k]["length"]):
        loci[k] = r
# the same gene recorded twice in one genome (two haplotype assemblies of one individual, a
# GenBank and a RefSeq annotation, or isoforms) is one row too: graded A-C representatives of
# one organism sharing >= 90 % of the shorter protein's 10-mers are merged into the longest.
# Not counted as independent support (no S2); D models are left apart (fusions share parts).
SEQ = {}
_name, _buf = None, []
_fa = os.path.join(D, "fasta", "unique_ge35992.fasta")
for _l in (open(_fa) if os.path.exists(_fa) else gzip.open(_fa + ".gz", "rt")):
    if _l[0] == ">":
        if _name:
            SEQ[_name] = "".join(_buf)
        _name, _buf = _l[1:].split()[0], []
    else:
        _buf.append(_l.strip())
if _name:
    SEQ[_name] = "".join(_buf)


def kmers(md5, k=10):
    q = SEQ.get(md5, "")
    return {q[i:i + k] for i in range(len(q) - k + 1)}


merged_into = {}
by_org = {}
for k, r in loci.items():
    if r["grade"] in "ABC" and r["md5"] in SEQ:
        by_org.setdefault(r["organism"], []).append(k)
for org, ks in by_org.items():
    if len(ks) < 2:
        continue
    ks.sort(key=lambda k: -int(loci[k]["length"]))
    km = {k: kmers(loci[k]["md5"]) for k in ks}
    for i, a in enumerate(ks):
        if a in merged_into:
            continue
        for b in ks[i + 1:]:
            if b in merged_into:
                continue
            shared = len(km[a] & km[b]) / max(1, min(len(km[a]), len(km[b])))
            if shared >= 0.90:
                merged_into[b] = (a, round(shared, 3))
same_record = {}
for b, (a, sh) in merged_into.items():
    same_record.setdefault(a, []).append(f"{loci[b]['rep_accession']} ({loci[b]['length']} aa, grade {loci[b]['grade']}, "
                                         f"{sh:.0%} shared 10-mers)")
    best[a] = min(best[a], best[b], key=ORDER.index)
rows = []
for k, r in loci.items():
    if k in merged_into:
        continue
    rows.append({"accession": r["rep_accession"], "organism": r["organism"], "superkingdom": r["superkingdom"],
                 "length_aa": int(r["length"]), "grade": r["grade"], "grade_stage1": r["grade_stage1"],
                 "best_grade_other_models_at_locus": best[k], "databases": r["databases"],
                 "evidence": evidence(r), "still_missing": missing(r), "locus_key": k,
                 "same_gene_other_assemblies": ";".join(a for a in copies_of[find(r["md5"])]
                                                       if a != r["rep_accession"]),
                 "same_gene_other_records": "; ".join(same_record.get(k, []))})
for m, r in s0.items():
    if r["triage"] == "T-excl":
        why = [x for x in ("synthetic/patent", "no-live-record", "partial") if x in r["flags"]]
        rows.append({"accession": r["rep_accession"], "organism": r["organism"], "superkingdom": r["superkingdom"],
                     "length_aa": int(r["length"]), "grade": "X", "grade_stage1": "X",
                     "best_grade_other_models_at_locus": "", "databases": r["databases"],
                     "evidence": "excluded: " + ", ".join(why), "still_missing": "", "locus_key": r["locus_key"],
                     "same_gene_other_assemblies": "", "same_gene_other_records": ""})
rows.sort(key=lambda x: -x["length_aa"])
rank = 0
for x in rows:
    if x["grade"] != "X":
        rank += 1
        x["rank"] = rank
    else:
        x["rank"] = ""
cols = ["rank", "accession", "organism", "superkingdom", "length_aa", "grade", "grade_stage1",
        "best_grade_other_models_at_locus", "same_gene_other_assemblies", "same_gene_other_records", "databases",
        "evidence", "still_missing", "locus_key"]
with open(os.path.join(D, "leaderboard.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for x in rows:
        f.write("\t".join(str(x[c]).replace("\t", " ") for c in cols) + "\n")
from collections import Counter  # noqa: E402
graded = [x for x in rows if x["grade"] != "X"]
print("loci", len(graded), Counter(x["grade"] for x in graded), "excluded rows", len(rows) - len(graded))
print("loci > 45,212 aa", Counter(x["grade"] for x in graded if x["length_aa"] > 45212))
print("loci stage1 -> stage2", Counter((x["grade_stage1"], x["grade"]) for x in graded))
for x in graded[:25]:
    print(x["rank"], x["accession"], x["organism"][:26], x["length_aa"], x["grade"], "|", x["evidence"][:110])
