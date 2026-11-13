"""Search database for the PXD060105 (Sycon ciliatum spicule matrix) re-search.

The authors' database (raw/Sponge_Sycon_2023_HBWS01_CDS.fasta, 99,773 transcript-derived proteins)
holds the 50,516-aa RefSeq model XP_065194412.1 only as partial fragments. Entries of which >= 50 % of
the 10-mers occur in XP_065194412.1 are treated as fragments of the same gene and replaced by the full
model (as QAT16263.1 was replaced by WP_128698900.1 for PXD025008); a reversed and a shuffled copy of
the model are added as controls. Trailing stop symbols are removed. Writes db/search_db.fasta and
db/dropped_fragments.tsv. Usage: python3 build_db.py
"""
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
K = 10


def read(p):
    out, name = [], None
    for line in open(p):
        if line.startswith(">"):
            name = line[1:].strip()
            out.append([name, []])
        else:
            out[-1][1].append(line.strip())
    return [(n, "".join(s).rstrip("*")) for n, s in out]


giant = read(os.path.join(HERE, "XP_065194412.1.fasta"))[0][1]
gk = {giant[i:i + K] for i in range(len(giant) - K + 1)}
kept, dropped = [], []
for n, s in read(os.path.join(HERE, "raw", "Sponge_Sycon_2023_HBWS01_CDS.fasta")):
    ks = [s[i:i + K] for i in range(len(s) - K + 1)]
    share = sum(1 for k in ks if k in gk) / len(ks) if ks else 0
    (dropped if share >= 0.5 else kept).append((n, s, share))
rng = random.Random(20261113)
shuf = list(giant)
rng.shuffle(shuf)
entries = [(n.split("|")[0] + " " + n, s) for n, s, _ in kept]
entries += [("XP_065194412.1 giant RefSeq model (Sycon ciliatum, 50,516 aa)", giant),
            ("CTRLREV_giant reversed XP_065194412.1", giant[::-1]),
            ("CTRLSHUF_giant shuffled XP_065194412.1 seed 20261113", "".join(shuf))]
os.makedirs(os.path.join(HERE, "db"), exist_ok=True)
with open(os.path.join(HERE, "db", "search_db.fasta"), "w") as f:
    for n, s in entries:
        f.write(f">{n}\n")
        for i in range(0, len(s), 80):
            f.write(s[i:i + 80] + "\n")
with open(os.path.join(HERE, "db", "dropped_fragments.tsv"), "w") as f:
    f.write("entry\tlength\tshare_of_10mers_in_XP_065194412.1\n")
    for n, s, sh in dropped:
        f.write(f"{n}\t{len(s)}\t{sh:.3f}\n")
print(f"kept {len(kept)}, dropped as fragments of the model {len(dropped)}, total entries {len(entries)}")
