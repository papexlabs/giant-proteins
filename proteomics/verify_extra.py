"""Independent re-computation of the numbers the proteomics fork reported beyond full6_summary.json.

Usage: python3 verify_extra.py
Reads raw/*.MGF, db/search_db.fasta and search/out/results.sage.tsv. Writes results/verify_extra.json:
  spectra in the six peak lists (BEGIN IONS count), database entries, best protein q of the giant and
  of its reversed and shuffled copies, theoretical fully tryptic peptides of 7-40 aa per 1,000-aa bin
  (K/R not before P, no missed cleavage), unique peptide starts in the first 6,000 aa, last residue
  covered, and internal repeat blocks of the giant (runs of 20-mers that occur more than once in it,
  merged when closer than 10 aa, blocks of >= 100 aa).
"""
import csv
import glob
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

n = 0
scans = {}
for f in sorted(glob.glob(os.path.join(HERE, "raw", "*.MGF"))):
    with open(f, "rb") as fh:
        for line in fh:
            if line.startswith(b"BEGIN IONS"):
                n += 1
            elif line.startswith(b"TITLE="):
                m = re.match(rb"TITLE=(.+)\.(\d+)\.\d+\.\d+\.dta", line.strip())
                scans.setdefault(m.group(1), set()).add(int(m.group(2)))
out["spectra_in_peak_lists"] = n
# share of scan numbers present: spectra per fraction / highest scan number in that fraction (scan numbers
# also count MS1 scans, so this is a lower bound on the share of MS2 spectra present)
out["fractions"] = len(scans)
out["scan_numbers_max_sum"] = sum(max(v) for v in scans.values())
out["share_of_scan_numbers_present"] = round(sum(len(v) for v in scans.values()) / out["scan_numbers_max_sum"], 3)

db, name = {}, None
for line in open(os.path.join(HERE, "db", "search_db.fasta")):
    if line.startswith(">"):
        name = line[1:].split()[0]
        db[name] = []
    else:
        db[name].append(line.strip())
db = {k: "".join(v) for k, v in db.items()}
out["database_entries"] = len(db)
GI = next(k for k in db if "WP_128698900.1" in k)
g = db[GI]

best = {"giant": 1.0, "CTRLREV_giant": 1.0, "CTRLSHUF_giant": 1.0}
with open(os.path.join(HERE, "search", "out", "results.sage.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["label"] != "1":
            continue
        prots = r["proteins"].split(";")
        for key, acc in (("giant", GI), ("CTRLREV_giant", "CTRLREV_giant"), ("CTRLSHUF_giant", "CTRLSHUF_giant")):
            if acc in prots:
                best[key] = min(best[key], float(r["protein_q"]))
out["best_protein_q"] = best

cuts = [0] + [m.end() for m in re.finditer(r"[KR](?!P)", g)] + [len(g)]
cuts = sorted(set(cuts))
bins = [0] * ((len(g) + 999) // 1000)
for a, b in zip(cuts, cuts[1:]):
    if 7 <= b - a <= 40:
        bins[a // 1000] += 1
out["theoretical_tryptic_7_40_per_1000aa"] = bins
out["theoretical_min_max_full_bins"] = (min(bins[:-1]), max(bins[:-1]))

pep = list(csv.DictReader(open(os.path.join(HERE, "results", "full6_giant_peptides.tsv")), delimiter="\t"))
uniq = [p for p in pep if p["unique_to_giant"] == "True"]
starts = [int(x) for p in uniq for x in str(p["start_positions"]).replace("[", "").replace("]", "").split(",") if x.strip()]
out["unique_peptides"] = len(uniq)
out["unique_starts_first_6000"] = sum(1 for s in starts if s <= 6000)
out["unique_starts_after_6000"] = sorted(s for s in starts if s > 6000)
ends = [s + len(p["peptide"]) - 1 for p in uniq for s in [int(x) for x in str(p["start_positions"]).replace("[", "").replace("]", "").split(",") if x.strip()]]
out["last_residue_covered"] = max(ends)
out["distal_peptides_psms"] = [(p["peptide"], int(p["start_positions"]), int(p["psms"]), float(p["peptide_q"]))
                               for p in uniq if int(str(p["start_positions"]).split(";")[0]) > 6000]
out["fraction_of_length_up_to_last_residue"] = round(max(ends) / len(g), 3)

K = 20
seen = {}
for i in range(len(g) - K + 1):
    seen.setdefault(g[i:i + K], []).append(i)
rep = sorted(i for v in seen.values() if len(v) > 1 for i in v)
blocks = []
for i in rep:
    if blocks and i <= blocks[-1][1] + 10:
        blocks[-1][1] = max(blocks[-1][1], i + K)
    else:
        blocks.append([i, i + K])
blocks = [(a + 1, b) for a, b in blocks if b - a >= 100]
out["repeat_blocks_ge100aa"] = blocks
out["repeat_span"] = (blocks[0][0], blocks[-1][1]) if blocks else None
out["unique_peptides_inside_or_after_first_block"] = sum(1 for s in starts if blocks and s >= blocks[0][0])

json.dump(out, open(os.path.join(HERE, "results", "verify_extra.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "theoretical_tryptic_7_40_per_1000aa"}, indent=1))
