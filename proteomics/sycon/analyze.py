"""Peptide evidence for XP_065194412.1 in the PXD060105 re-search (search/out, Sage).

Thresholds as for PXD025008: spectrum, peptide and protein q <= 0.01. 'Unique' = the peptide (I = L) occurs in
no other forward entry of db/search_db.fasta. Comparison proteins, fixed from the authors' own report before
our results were seen: their top identification (CAH0794662, Calcarin 3, most exclusive unique peptides) and
two housekeeping proteins in their list (CAH0809064 elongation factor 2-like, CAH0788721 tubulin beta-like);
plus the share of the authors' 38 reported proteins that this search also finds at 1 % protein FDR.
Writes results/summary.json and results/giant_peptides.tsv. Usage: python3 analyze.py
"""
import csv
import json
import os
import re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
Q = 0.01
IL = str.maketrans("I", "L")
db, name = {}, None
for line in open(os.path.join(HERE, "db", "search_db.fasta")):
    if line.startswith(">"):
        name = line[1:].split()[0]
        db[name] = []
    else:
        db[name].append(line.strip())
db = {k: "".join(v) for k, v in db.items()}
G = "XP_065194412.1"
giant = db[G]
others = "#".join(v.translate(IL) for k, v in db.items() if k != G and not k.startswith("CTRL"))
strip = lambda p: re.sub(r"\[.*?\]|[^A-Z]", "", p)  # noqa: E731
tgt, dec = [], []
with open(os.path.join(HERE, "search", "out", "results.sage.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        (tgt if r["label"] == "1" else dec).append(r)
ok = [r for r in tgt if float(r["spectrum_q"]) <= Q and float(r["peptide_q"]) <= Q and float(r["protein_q"]) <= Q]
prot_q = defaultdict(lambda: 1.0)
for r in tgt:
    for p in r["proteins"].split(";"):
        prot_q[p] = min(prot_q[p], float(r["protein_q"]))
summ = {"spectra_searched_ms2": None, "psms_reported": len(tgt) + len(dec),
        "target_psms_q01": sum(1 for r in tgt if float(r["spectrum_q"]) <= Q),
        "target_peptides_q01": len({strip(r["peptide"]) for r in tgt if float(r["peptide_q"]) <= Q}),
        "target_proteins_q01": len({p for p, q in prot_q.items() if q <= Q})}
md = [json.load(open(os.path.join(HERE, "mzml", f))) for f in os.listdir(os.path.join(HERE, "mzml")) if f.endswith("-metadata.json")]
summ["spectra_searched_ms2"] = sum(int(x["value"]) for d in md for x in d["MsData"] if x["name"] == "Number of MS2 spectra")


def entry_stats(acc):
    seq = db[acc]
    peps = defaultdict(lambda: {"psms": 0, "q": 1.0})
    for r in ok:
        if acc in r["proteins"].split(";"):
            p = strip(r["peptide"])
            peps[p]["psms"] += 1
            peps[p]["q"] = min(peps[p]["q"], float(r["peptide_q"]))
    s_il = seq.translate(IL)
    others_x = "#".join(v.translate(IL) for k, v in db.items() if k != acc and not k.startswith("CTRL"))
    rows = []
    for p, v in peps.items():
        pil = p.translate(IL)
        pos = [m.start() + 1 for m in re.finditer(f"(?={re.escape(pil)})", s_il)]
        rows.append({"peptide": p, "start_positions": ";".join(map(str, pos)), "psms": v["psms"],
                     "peptide_q": v["q"], "unique": pil not in others_x})
    cov = set()
    for x in rows:
        if x["unique"]:
            for a in x["start_positions"].split(";"):
                if a:
                    cov.update(range(int(a), int(a) + len(x["peptide"])))
    return rows, {"length": len(seq), "peptides_q01": len(rows), "unique_peptides_q01": sum(1 for x in rows if x["unique"]),
                  "psms_q01": sum(x["psms"] for x in rows), "coverage_by_unique_frac": round(len(cov) / len(seq), 3),
                  "protein_q": prot_q.get(acc, 1.0)}


grows, gst = entry_stats(G)
uniq = [x for x in grows if x["unique"]]
starts = sorted(int(a) for x in uniq for a in x["start_positions"].split(";") if a)
gst["fifths_with_unique_peptides"] = len({min(4, (s - 1) * 5 // len(giant)) for s in starts})
gst["first_unique_start"] = starts[0] if starts else None
gst["last_unique_start"] = starts[-1] if starts else None
summ["giant"] = gst
summ["controls_giant_copies"] = {k: {"psms_q01": sum(1 for r in tgt if k in r["proteins"].split(";") and float(r["spectrum_q"]) <= Q),
                                     "best_protein_q": prot_q.get(k, 1.0)}
                                 for k in ("CTRLREV_giant", "CTRLSHUF_giant")}
comp = {}
for acc in ("CAH0794662", "CAH0809064", "CAH0788721"):
    full = next(k for k in db if k == acc)
    comp[acc] = entry_stats(full)[1]
summ["comparison_proteins"] = comp
L = open(os.path.join(HERE, "raw", "Protein_Report_for_Sycon_newDB_mergeP_April_1123.xls")).read().splitlines()
i = next(k for k, line in enumerate(L) if line.startswith("Experiment name"))
h = L[i].split("\t")
rep = [dict(zip(h, line.split("\t"))) for line in L[i + 1:] if line.startswith("Sycon_newDB")]
found = 0
for r in rep:
    accs = [a.strip() for a in r["Protein accession numbers"].split(",")]
    if any(prot_q.get(a, 1.0) <= Q for a in accs):
        found += 1
summ["authors_reported_proteins"] = len(rep)
summ["database_entries_authors"] = sum(1 for line in open(os.path.join(HERE, "raw", "Sponge_Sycon_2023_HBWS01_CDS.fasta"))
                                       if line.startswith(">"))
summ["database_entries_searched"] = len(db)
summ["raw_files"] = len([f for f in os.listdir(os.path.join(HERE, "mzml")) if f.endswith(".mzML")])
summ["authors_proteins_found_q01"] = found
# This dataset is small (the authors report 37 proteins, most with 1-5 peptides): with ~60 confident
# peptides one decoy already puts the peptide-level q above 0.01, so the triple 1 % filter can pass
# nothing. Reported in addition, and labelled as such: spectrum-level 1 % FDR (PSM q <= 0.01) only.
sp = [r for r in tgt if float(r["spectrum_q"]) <= Q]
summ["min_peptide_q_any_target"] = min(float(r["peptide_q"]) for r in tgt)
summ["min_protein_q_any_target"] = min(float(r["protein_q"]) for r in tgt)


def sp_stats(acc):
    rs = [r for r in sp if acc in r["proteins"].split(";")]
    return {"psms_spectrum_q01": len(rs), "peptides_spectrum_q01": len({strip(r["peptide"]) for r in rs}),
            "psms_any_q": sum(1 for r in tgt if acc in r["proteins"].split(";"))}


summ["spectrum_level"] = {"target_psms": len(sp), "target_peptides": len({strip(r["peptide"]) for r in sp}),
                          "giant": sp_stats(G), "CTRLREV_giant": sp_stats("CTRLREV_giant"),
                          "CTRLSHUF_giant": sp_stats("CTRLSHUF_giant"),
                          "comparison_proteins": {a: sp_stats(a) for a in ("CAH0794662", "CAH0809064", "CAH0788721")}}
f_sp = 0
for r in rep:
    accs = [a.strip() for a in r["Protein accession numbers"].split(",")]
    if any(sp_stats(a)["psms_spectrum_q01"] > 0 for a in accs if a in db):
        f_sp += 1
summ["spectrum_level"]["authors_proteins_with_psm"] = f_sp
with open(os.path.join(HERE, "results", "giant_peptides.tsv"), "w") as f:
    f.write("peptide\tstart_positions\tpsms\tpeptide_q\tunique_to_giant\n")
    for x in sorted(grows, key=lambda x: int(x["start_positions"].split(";")[0] or 0)):
        f.write(f"{x['peptide']}\t{x['start_positions']}\t{x['psms']}\t{x['peptide_q']}\t{x['unique']}\n")
json.dump(summ, open(os.path.join(HERE, "results", "summary.json"), "w"), indent=1)
print(json.dumps(summ, indent=1))
