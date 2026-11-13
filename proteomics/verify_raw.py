"""Checks of the independent RAW-file search of PXD025008 (search/out_raw), and its comparison
with the deposited-peak-list search (search/out).

Usage: python3 verify_raw.py
Reads mzml/*-metadata.json (ThermoRawFileParser), mzml/convert.log, search/out_raw/results.sage.tsv,
results/raw60_giant_peptides.tsv, results/full6_giant_peptides.tsv, results/verify_extra.json and the six
MGF peak lists. Writes results/verify_raw.json.
"""
import csv
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

ms2 = ms1 = 0
inst = set()
for f in sorted(glob.glob(os.path.join(HERE, "mzml", "*-metadata.json"))):
    d = json.load(open(f))
    for x in d.get("MsData", []):
        if x["name"] == "Number of MS2 spectra":
            ms2 += int(x["value"])
        if x["name"] == "Number of MS1 spectra":
            ms1 += int(x["value"])
    for x in d.get("InstrumentProperties", []):
        if x["name"] == "Thermo Scientific instrument model":
            inst.add(x["value"])
out["raw_files"] = len(glob.glob(os.path.join(HERE, "mzml", "*.mzML")))
out["ms2_spectra_raw"] = ms2
out["ms1_spectra_raw"] = ms1
out["instrument"] = sorted(inst)
log = open(os.path.join(HERE, "mzml", "convert.log")).read()
out["conversions_completed"] = len(re.findall(r"Processing completed (\d+) errors", log))
out["conversion_errors"] = sum(int(x) for x in re.findall(r"Processing completed (\d+) errors", log))
mgf = 0
for f in glob.glob(os.path.join(HERE, "raw", "*.MGF")):
    with open(f, "rb") as fh:
        mgf += sum(1 for line in fh if line.startswith(b"BEGIN IONS"))
out["mgf_spectra"] = mgf
out["mgf_share_of_raw_ms2"] = round(mgf / ms2, 3)

GI = "Velamenicoccus|WP_128698900.1|BU251_00165"
best = {"giant": 1.0, "CTRLREV_giant": 1.0, "CTRLSHUF_giant": 1.0}
with open(os.path.join(HERE, "search", "out_raw", "results.sage.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["label"] != "1":
            continue
        prots = r["proteins"].split(";")
        for key, acc in (("giant", GI), ("CTRLREV_giant", "CTRLREV_giant"), ("CTRLSHUF_giant", "CTRLSHUF_giant")):
            if acc in prots:
                best[key] = min(best[key], float(r["protein_q"]))
out["best_protein_q"] = best


def peps(name):
    rows = list(csv.DictReader(open(os.path.join(HERE, "results", name)), delimiter="\t"))
    return [r for r in rows if r["unique_to_giant"] == "True"]


raw, mgfp = peps("raw60_giant_peptides.tsv"), peps("full6_giant_peptides.tsv")
starts = [(int(p["start_positions"].split(";")[0]), p) for p in raw]
out["raw_unique_peptides"] = len(raw)
out["raw_starts_first_6000"] = sum(1 for s, _ in starts if s <= 6000)
out["raw_after_6000"] = sorted((s, int(p["psms"]), float(p["peptide_q"]), p["peptide"]) for s, p in starts if s > 6000)
out["raw_last_residue"] = max(s + len(p["peptide"]) - 1 for s, p in starts)
L = 39686
out["raw_last_residue_pct"] = round(100 * out["raw_last_residue"] / L, 1)
a, b = {p["peptide"] for p in raw}, {p["peptide"] for p in mgfp}
out["peptides_in_both_searches"] = len(a & b)
out["peptides_only_raw"] = len(a - b)
out["peptides_only_mgf"] = len(b - a)
out["mgf_only_after_6000"] = sorted(int(p["start_positions"].split(";")[0]) for p in mgfp
                                    if p["peptide"] in b - a and int(p["start_positions"].split(";")[0]) > 6000)
vx = json.load(open(os.path.join(HERE, "results", "verify_extra.json")))
blocks = vx["repeat_blocks_ge100aa"]
out["raw_peptides_inside_or_after_repeats"] = sum(1 for s, _ in starts if s >= blocks[0][0])
json.dump(out, open(os.path.join(HERE, "results", "verify_raw.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
