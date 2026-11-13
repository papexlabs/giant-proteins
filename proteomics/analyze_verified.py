"""Peptide evidence for the giant WP_128698900.1 in the PXD025008 reanalysis (Sage output).

Usage: python3 analyze_verified.py <sage output dir> <label>
Reads <dir>/results.sage.tsv and db/search_db.fasta. Thresholds: PSM, peptide and protein
q-values from Sage (target-decoy) all <= 0.01. 'Unique' = the peptide (I = L) occurs in no other
forward entry of the search database (checked here on the sequences, not on Sage's protein list).
Writes results/<label>_summary.json, results/<label>_giant_peptides.tsv, results/<label>_controls.tsv
"""
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR, LABEL = sys.argv[1], sys.argv[2]
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)
Q = 0.01
IL = str.maketrans("I", "L")

db, name = {}, None
for line in open(os.path.join(HERE, "db", "search_db.fasta")):
    if line.startswith(">"):
        name = line[1:].strip()
        db[name] = []
    else:
        db[name].append(line.strip())
db = {k: "".join(v) for k, v in db.items()}
first = {k.split()[0]: k for k in db}
GIANT = next(k for k in db if "WP_128698900.1" in k)
giant = db[GIANT]
others = "#".join(v.translate(IL) for k, v in db.items() if k != GIANT and not k.startswith("CTRL"))
controls = {"reversed giant": next(k for k in db if k.startswith("CTRLREV")),
            "shuffled giant": next(k for k in db if k.startswith("CTRLSHUF"))}
# abundant normal proteins of the same organism, chosen by name before looking at the results
# (the search database keeps only 'Velamenicoccus|<protein_id>|<locus_tag>'; product names come from
# the GenBank CDS file of the closed genome, db/CP019384.1_cds_aa.fasta)
NORMAL = {"translation elongation factor Tu": "elongation factor Tu", "50S ribosomal protein L7/L12": "ribosomal protein L7/L12",
          "chaperonin GroEL": "chaperonin GroEL"}
pid_of = {}
for line in open(os.path.join(HERE, "db", "CP019384.1_cds_aa.fasta")):
    if line.startswith(">") and "[protein=" in line and "[protein_id=" in line:
        prod = line.split("[protein=")[1].split("]")[0]
        pid_of.setdefault(prod, line.split("[protein_id=")[1].split("]")[0])
normal = {}
for prod, lab in NORMAL.items():
    if prod in pid_of:
        hits = [k for k in db if k.startswith(f"Velamenicoccus|{pid_of[prod]}|")]
        if hits:
            normal[lab] = hits[0]


def strip(p):
    return re.sub(r"\[[^\]]*\]|[^A-Z]", "", p)


rows = list(csv.DictReader(open(os.path.join(HERE, OUTDIR, "results.sage.tsv")), delimiter="\t"))
tgt = [r for r in rows if r["label"] == "1"]
dec = [r for r in rows if r["label"] == "-1"]
ok = [r for r in tgt if float(r["spectrum_q"]) <= Q and float(r["peptide_q"]) <= Q and float(r["protein_q"]) <= Q]
summ = {"sage_output": OUTDIR, "files": sorted({r["filename"] for r in rows}), "psms_reported": len(rows),
        "target_psms_q01": sum(1 for r in tgt if float(r["spectrum_q"]) <= Q),
        "decoy_psms_q01": sum(1 for r in dec if float(r["spectrum_q"]) <= Q),
        "target_peptides_q01": len({strip(r["peptide"]) for r in tgt if float(r["peptide_q"]) <= Q}),
        "target_proteins_q01": len({p for r in tgt if float(r["protein_q"]) <= Q for p in r["proteins"].split(";")})}


def protein_hits(entry):
    seq = db[entry]
    rest = others if entry == GIANT else "#".join(v.translate(IL) for k, v in db.items()
                                                  if k != entry and not k.startswith("CTRL"))
    peps = defaultdict(lambda: {"psms": 0, "q": 1.0})
    for r in ok:
        prots = [first.get(p, p) for p in r["proteins"].split(";")]
        if entry not in prots:
            continue
        pep = strip(r["peptide"])
        peps[pep]["psms"] += 1
        peps[pep]["q"] = min(peps[pep]["q"], float(r["peptide_q"]))
    out = []
    s_il = seq.translate(IL)
    for pep, v in peps.items():
        p_il = pep.translate(IL)
        pos = [m.start() + 1 for m in re.finditer(f"(?={re.escape(p_il)})", s_il)]
        out.append({"peptide": pep, "psms": v["psms"], "peptide_q": v["q"], "positions": pos,
                    "unique_to_entry": p_il not in rest})
    return sorted(out, key=lambda x: (x["positions"] or [0])[0])


gp = protein_hits(GIANT)
uniq = [p for p in gp if p["unique_to_entry"]]
cov = bytearray(len(giant))
for p in uniq:
    for a in p["positions"]:
        cov[a - 1:a - 1 + len(p["peptide"])] = b"\x01" * len(p["peptide"])
regions, i = [], 0
while i < len(cov):
    if cov[i]:
        j = i
        while j + 1 < len(cov) and cov[j + 1]:
            j += 1
        regions.append((i + 1, j + 1))
        i = j + 1
    else:
        i += 1
starts = sorted(a for p in uniq for a in p["positions"])
summ["giant"] = {"entry": GIANT, "length": len(giant), "peptides_q01": len(gp), "unique_peptides_q01": len(uniq),
                 "shared_peptides_q01": len(gp) - len(uniq),
                 "unique_psms_q01": sum(p["psms"] for p in uniq),
                 "unique_peptides_at_several_positions_in_the_giant": sum(1 for p in uniq if len(p["positions"]) > 1),
                 "residues_covered_by_unique_peptides": sum(cov), "coverage_frac": round(sum(cov) / len(giant), 4),
                 "covered_regions": len(regions), "first_unique_position": starts[0] if starts else None,
                 "last_unique_position": starts[-1] if starts else None}
fifths = Counter(min(4, (a - 1) * 5 // len(giant)) for a in starts)
summ["giant"]["unique_peptides_per_fifth_of_length"] = [fifths.get(k, 0) for k in range(5)]
summ["controls_giant_copies"] = {}
for lab, k in controls.items():
    anyq = [r for r in tgt if k in [first.get(p, p) for p in r["proteins"].split(";")]]
    summ["controls_giant_copies"][lab] = {"psms_q01": sum(1 for r in ok if k in [first.get(p, p) for p in r["proteins"].split(";")]),
                                          "psms_any_q": len(anyq)}
summ["normal_proteins"] = {}
for lab, k in normal.items():
    ph = protein_hits(k)
    c = bytearray(len(db[k]))
    for p in ph:
        for a in p["positions"]:
            c[a - 1:a - 1 + len(p["peptide"])] = b"\x01" * len(p["peptide"])
    cu = bytearray(len(db[k]))
    for p in ph:
        if p["unique_to_entry"]:
            for a in p["positions"]:
                cu[a - 1:a - 1 + len(p["peptide"])] = b"\x01" * len(p["peptide"])
    summ["normal_proteins"][lab] = {"entry": k.split(" [")[0], "length": len(db[k]), "peptides_q01": len(ph),
                                    "unique_peptides_q01": sum(1 for p in ph if p["unique_to_entry"]),
                                    "psms_q01": sum(p["psms"] for p in ph), "coverage_frac": round(sum(c) / len(db[k]), 3),
                                    "coverage_by_unique_frac": round(sum(cu) / len(db[k]), 3)}
json.dump(summ, open(os.path.join(RES, f"{LABEL}_summary.json"), "w"), indent=1)
with open(os.path.join(RES, f"{LABEL}_giant_peptides.tsv"), "w") as f:
    f.write("peptide\tstart_positions\tlength\tpsms\tpeptide_q\tunique_to_giant\n")
    for p in gp:
        f.write(f"{p['peptide']}\t{';'.join(map(str, p['positions']))}\t{len(p['peptide'])}\t{p['psms']}\t{p['peptide_q']:.4g}\t{p['unique_to_entry']}\n")
print(json.dumps(summ, indent=1))
