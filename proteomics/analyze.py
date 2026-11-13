"""Summarise the Sage search of PXD025008 for the giant WP_128698900.1 and the controls.

Reads search/out/results.sage.tsv, db/search_db.fasta and the peak lists in raw/; writes
results/summary.json, results/giant_peptides.tsv, results/coverage_map.tsv and results/control_proteins.tsv.
Peptide-level and protein-level q-values are Sage's (picked-peptide and picked-protein target-decoy).
"Unique to the giant" = the peptide sequence (I and L treated as equal) occurs in no other protein of the
search database, checked here independently of Sage's protein list.
"""
import csv
import glob
import json
import os
import re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)
GIANT = "Velamenicoccus|WP_128698900.1|BU251_00165"
Q = 0.01


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    seqs[name] = "".join(buf)
    return seqs


db = read_fasta(os.path.join(HERE, "db", "search_db.fasta"))
key = {n.split()[0]: n for n in db}                 # Sage reports the first word of the header
gname = next(n for n in db if n.startswith(GIANT))
giant = db[gname]
IL = str.maketrans("I", "L")
allseq_other = "#".join(s.translate(IL) for n, s in db.items() if n != gname)
giant_il = giant.translate(IL)

rows = list(csv.DictReader(open(os.path.join(HERE, "search", "out", "results.sage.tsv")), delimiter="\t"))
target = [r for r in rows if r["label"] == "1"]
n_spectra = 0
for p in sorted(glob.glob(os.path.join(HERE, "raw", "Mudpit_*.MGF"))):
    with open(p, errors="replace") as f:
        n_spectra += sum(1 for line in f if line.startswith("BEGIN IONS"))
psm1 = [r for r in target if float(r["spectrum_q"]) <= Q]
pep1 = {re.sub(r"\[.*?\]", "", r["peptide"]) for r in target if float(r["peptide_q"]) <= Q}
prot1 = set()
for r in target:
    if float(r["protein_q"]) <= Q:
        prot1.update(r["proteins"].split(";"))


def strip(p):
    return re.sub(r"\[.*?\]|[^A-Z]", "", p)


def protein_hits(acc_prefix):
    """PSMs (target, peptide q <= 1 %) whose protein list contains a protein starting with acc_prefix."""
    hits = defaultdict(lambda: {"psms": 0, "best_peptide_q": 1.0, "best_spectrum_q": 1.0, "protein_q": 1.0,
                                "num_proteins": 0, "proteins": ""})
    for r in target:
        prots = r["proteins"].split(";")
        if not any(p.startswith(acc_prefix) for p in prots):
            continue
        if float(r["peptide_q"]) > Q:
            continue
        h = hits[strip(r["peptide"])]
        h["psms"] += 1
        h["best_peptide_q"] = min(h["best_peptide_q"], float(r["peptide_q"]))
        h["best_spectrum_q"] = min(h["best_spectrum_q"], float(r["spectrum_q"]))
        h["protein_q"] = min(h["protein_q"], float(r["protein_q"]))
        h["num_proteins"] = int(r["num_proteins"])
        h["proteins"] = r["proteins"]
    return hits


def positions(seq_il, pep):
    p = pep.translate(IL)
    out, i = [], seq_il.find(p)
    while i >= 0:
        out.append(i + 1)
        i = seq_il.find(p, i + 1)
    return out


# ---- the giant ---------------------------------------------------------------------------------
gh = protein_hits("Velamenicoccus|WP_128698900.1")
gtab = []
for pep, h in sorted(gh.items(), key=lambda x: positions(giant_il, x[0])[0] if positions(giant_il, x[0]) else 0):
    pos = positions(giant_il, pep)
    unique = pep.translate(IL) not in allseq_other
    gtab.append({"peptide": pep, "start": pos[0] if pos else "", "end": pos[0] + len(pep) - 1 if pos else "",
                 "occurrences_in_giant": len(pos), "all_starts": ";".join(map(str, pos)), "psms_q1": h["psms"],
                 "best_peptide_q": h["best_peptide_q"], "protein_q": h["protein_q"],
                 "unique_to_giant": "yes" if unique else "no", "sage_num_proteins": h["num_proteins"]})
with open(os.path.join(OUT, "giant_peptides.tsv"), "w") as f:
    cols = list(gtab[0].keys()) if gtab else ["peptide"]
    f.write("\t".join(cols) + "\n")
    for t in gtab:
        f.write("\t".join(str(t[c]) for c in cols) + "\n")
uniq = [t for t in gtab if t["unique_to_giant"] == "yes"]
# coverage map of unique peptides (all occurrences inside the giant)
iv = sorted((s, s + len(t["peptide"]) - 1) for t in uniq for s in map(int, t["all_starts"].split(";")))
merged = []
for a, b in iv:
    if merged and a <= merged[-1][1] + 1:
        merged[-1][1] = max(merged[-1][1], b)
    else:
        merged.append([a, b])
covered = sum(b - a + 1 for a, b in merged)
with open(os.path.join(OUT, "coverage_map.tsv"), "w") as f:
    f.write("start\tend\tlength\n")
    for a, b in merged:
        f.write(f"{a}\t{b}\t{b - a + 1}\n")
bins = [0] * ((len(giant) + 999) // 1000)
for t in uniq:
    for s in map(int, t["all_starts"].split(";")):
        bins[(s - 1) // 1000] += 1
# internal repeats of the giant: stretches whose 20-mers occur more than once in the protein
K = 20
cnt = defaultdict(list)
for i in range(len(giant) - K + 1):
    cnt[giant_il[i:i + K]].append(i)
rep = [0] * len(giant)
for v in cnt.values():
    if len(v) > 1:
        for i in v:
            for j in range(i, i + K):
                rep[j] = 1
blocks, i = [], 0
while i < len(rep):
    if rep[i]:
        j = i
        while j + 1 < len(rep) and rep[j + 1]:
            j += 1
        if j - i + 1 >= 100:
            blocks.append([i + 1, j + 1])
        i = j + 1
    else:
        i += 1
starts_u = sorted(int(s) for t in uniq for s in t["all_starts"].split(";") if t["occurrences_in_giant"] == 1)
around = []
for a, b in blocks:
    around.append({"repeat_block": [a, b], "unique_single_copy_peptides_before": sum(1 for s in starts_u if s < a),
                   "unique_single_copy_peptides_after": sum(1 for s in starts_u if s > b),
                   "peptides_inside": sum(1 for t in uniq for s in map(int, t["all_starts"].split(";")) if a <= s <= b)})

# ---- controls ---------------------------------------------------------------------------------
ctrl = {}
for tag in ("CTRLREV_giant", "CTRLSHUF_giant"):
    h = protein_hits(tag)
    psm_any = sum(1 for r in target if tag in r["proteins"] and float(r["spectrum_q"]) <= Q)
    ctrl[tag] = {"peptides_q1": len(h), "psms_spectrum_q1": psm_any,
                 "protein_q": min([float(r["protein_q"]) for r in target if tag in r["proteins"]], default=None)}
# abundant normal proteins of the organism: named housekeeping proteins + the most-detected proteins
names = {n.split()[0]: n for n in db}
genome_hdr = {}
for line in open(os.path.join(HERE, "db", "CP019384.1_cds_aa.fasta")):
    if line.startswith(">"):
        pid = line.split("[protein_id=")[1].split("]")[0] if "[protein_id=" in line else ""
        prot = line.split("[protein=")[1].split("]")[0] if "[protein=" in line else ""
        genome_hdr[pid] = prot
psm_by_prot = defaultdict(int)
for r in psm1:
    for p in r["proteins"].split(";"):
        if p.startswith("Velamenicoccus|"):
            psm_by_prot[p] += 1
top = sorted(psm_by_prot.items(), key=lambda x: -x[1])[:12]
want = [p for p, _ in top]
for pid, prot in genome_hdr.items():
    if re.search(r"elongation factor Tu|chaperonin GroEL|molecular chaperone GroEL", prot, re.I):
        want.append(next((n.split()[0] for n in db if n.startswith(f"Velamenicoccus|{pid}|")), ""))
ctab = []
for acc in dict.fromkeys(w for w in want if w):
    full = key.get(acc, acc)
    seq = db[full]
    s_il = seq.translate(IL)
    h = protein_hits(acc)
    up = [p for p in h if p.translate(IL) not in "#".join(x.translate(IL) for n, x in db.items() if n != full)]
    cov = [0] * len(seq)
    for p in up:
        for s in positions(s_il, p):
            for k in range(s - 1, s - 1 + len(p)):
                cov[k] = 1
    pid = acc.split("|")[1]
    ctab.append({"protein": acc, "name": "giant" if pid == "WP_128698900.1" else genome_hdr.get(pid, ""),
                 "length": len(seq), "psms_q1": psm_by_prot.get(acc, 0), "peptides_q1": len(h),
                 "unique_peptides_q1": len(up), "unique_coverage_pct": round(100 * sum(cov) / len(seq), 1),
                 "protein_q": min([v["protein_q"] for v in h.values()], default=None)})
with open(os.path.join(OUT, "control_proteins.tsv"), "w") as f:
    cols = list(ctab[0].keys())
    f.write("\t".join(cols) + "\n")
    for t in ctab:
        f.write("\t".join(str(t[c]) for c in cols) + "\n")

summary = {
    "dataset": "PRIDE PXD025008 (Kizina et al. 2022, Appl Environ Microbiol 88:e02407-21)",
    "inputs": {"peak_lists": sorted(os.path.basename(p) for p in glob.glob(os.path.join(HERE, "raw", "Mudpit_*.MGF"))),
               "database": "db/search_db.fasta (built by build_db.py)", "search": "run_search.sh (Sage)",
               "config": "search/sage_config.json"},
    "spectra_in_peak_lists": n_spectra,
    "target_psms_spectrum_q_le_0.01": len(psm1),
    "target_peptides_peptide_q_le_0.01": len(pep1),
    "target_proteins_protein_q_le_0.01": len(prot1),
    "giant": {"accession": "WP_128698900.1", "length": len(giant),
              "protein_q": min([t["protein_q"] for t in gtab], default=None),
              "peptides_q1_any": len(gtab), "peptides_q1_unique_to_giant": len(uniq),
              "psms_on_unique_peptides": sum(t["psms_q1"] for t in uniq),
              "unique_peptides_occurring_more_than_once_in_giant": sum(1 for t in uniq if t["occurrences_in_giant"] > 1),
              "covered_residues_by_unique_peptides": covered,
              "coverage_pct": round(100 * covered / len(giant), 2),
              "distinct_regions": len(merged),
              "first_unique_peptide_start": iv[0][0] if iv else None,
              "last_unique_peptide_end": iv[-1][1] if iv else None,
              "unique_peptide_starts_per_1000aa_bin": bins,
              "internal_repeat_blocks_20mer": blocks,
              "peptides_around_repeat_blocks": around},
    "controls": ctrl,
    "abundant_proteins_same_organism": ctab,
}
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != "abundant_proteins_same_organism"}, indent=1)[:4000])
for t in ctab:
    print(t)
