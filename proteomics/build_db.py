"""Search database for the PXD025008 reanalysis.

Inputs (downloaded by run_download.sh):
  raw/OP3_AA_20170330_fw_rev_cont.fasta   the authors' metaproteome database (forward + REVERSE_ + contaminants)
  db/CP019384.1_cds_aa.fasta              closed genome of Candidatus Velamenicoccus archaeovorus (GenBank CDS)
  db/WP_128698900.1.fasta                 the giant, RefSeq version (39,686 aa; GenBank QAT16263.1 = same
                                           sequence without the first 9 residues)
Output: db/search_db.fasta
  - every forward entry of the authors' database except the 2017 'Omnitrophica_OP3_*' bin, which is
    replaced by the closed-genome proteome (co-cultured community and contaminants kept as they are;
    the authors' REVERSE_ decoys are dropped because Sage makes its own)
  - closed-genome proteome with QAT16263.1 replaced by WP_128698900.1
  - two controls: CTRLREV_giant (the giant reversed) and CTRLSHUF_giant (the giant shuffled, seed 20261113)
"""
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))


def read(path):
    out, name = [], None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
            out.append([name, []])
        elif line.strip():
            out[-1][1].append(line.strip())
    return [(n, "".join(s)) for n, s in out]


authors = read(os.path.join(HERE, "raw", "OP3_AA_20170330_fw_rev_cont.fasta"))
genome = read(os.path.join(HERE, "db", "CP019384.1_cds_aa.fasta"))
giant = read(os.path.join(HERE, "db", "WP_128698900.1.fasta"))[0][1]

keep = [(n, s) for n, s in authors if not n.startswith("REVERSE_") and not n.startswith("Omnitrophica_OP3")]
dropped_op3 = sum(1 for n, _ in authors if n.startswith("Omnitrophica_OP3"))
genome_entries = []
for n, s in genome:
    pid = n.split("[protein_id=")[1].split("]")[0] if "[protein_id=" in n else n.split()[0]
    lt = n.split("[locus_tag=")[1].split("]")[0] if "[locus_tag=" in n else ""
    if pid == "QAT16263.1":
        continue
    genome_entries.append((f"Velamenicoccus|{pid}|{lt}", s))
genome_entries.append(("Velamenicoccus|WP_128698900.1|BU251_00165 giant (RefSeq version of QAT16263.1)", giant))
rnd = random.Random(20261113)
shuf = list(giant)
rnd.shuffle(shuf)
controls = [("CTRLREV_giant reversed WP_128698900.1", giant[::-1]),
            ("CTRLSHUF_giant shuffled WP_128698900.1 seed 20261113", "".join(shuf))]
allp = keep + genome_entries + controls
with open(os.path.join(HERE, "db", "search_db.fasta"), "w") as f:
    for n, s in allp:
        f.write(f">{n}\n")
        for i in range(0, len(s), 80):
            f.write(s[i:i + 80] + "\n")
print(f"authors' forward entries kept {len(keep)} (dropped 2017 OP3 bin: {dropped_op3}); "
      f"closed-genome proteins {len(genome_entries)}; controls {len(controls)}; total {len(allp)}; "
      f"residues {sum(len(s) for _, s in allp)}")
