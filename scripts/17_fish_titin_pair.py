"""Follow-up (a): are the merged fish 'titin' models exactly zebrafish's tandem ttna + ttnb pair?

Each 1,000-aa window of the two merged models (Cirrhinus KAK2904762.1, Bagarius TSR51413.1)
is searched with phmmer against zebrafish ttna (XP_087076962.1, 31,787 aa) and ttnb
(XP_073769224.1, 26,697 aa); per window the better-scoring paralogue is recorded. Zebrafish
order (NCBI Gene, GRCz12tu chr9 NC_133184.1, both genes on the minus strand): ttna at
44.63-44.81 Mb lies upstream (5') of ttnb at 44.47-44.61 Mb, 21 kb apart, so a read-through
model of the pair would be ttna followed by ttnb.
Output: data/followups/fish_titin_pair.tsv and a printed summary.
"""
import csv
import os
import pyhmmer

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
alpha = pyhmmer.easel.Alphabet.amino()


def fasta(path):
    out, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


zf = fasta(os.path.join(D, "followups", "zebrafish_ttna_ttnb.fasta"))
label = {"XP_087076962.1": "ttna", "XP_073769224.1": "ttnb"}
targets = pyhmmer.easel.DigitalSequenceBlock(alpha, [
    pyhmmer.easel.TextSequence(name=k.encode(), sequence=v).digitize(alpha) for k, v in zf.items()])
checks = {r["rep_accession"]: r for r in csv.DictReader(open(os.path.join(D, "locus_checks.tsv")), delimiter="\t")}
uniq = {}
name = None
for line in open(os.path.join(D, "fasta", "unique_ge35992.fasta")):
    if line.startswith(">"):
        name = line[1:].split()[0]
        uniq[name] = []
    else:
        uniq[name].append(line.strip())
rows = []
for acc in ("KAK2904762.1", "TSR51413.1"):
    seq = "".join(uniq[checks[acc]["md5"]])
    W = 1000
    qs = [pyhmmer.easel.TextSequence(name=f"{acc}|{o}".encode(), sequence=seq[o:o + W]).digitize(alpha)
          for o in range(0, len(seq), W) if len(seq[o:o + W]) >= 200]
    hits = {}
    for th in pyhmmer.hmmer.phmmer(qs, targets, cpus=6, E=1e-5):
        q = th.query.name.decode() if isinstance(th.query.name, bytes) else th.query.name
        for h in th:
            t = h.name.decode() if isinstance(h.name, bytes) else h.name
            dom = max(h.domains, key=lambda d: d.score)
            al = dom.alignment
            hits.setdefault(q, {})[label[t]] = (round(h.score, 1), al.target_from, al.target_to)
    for q in sorted(hits, key=lambda x: int(x.split("|")[1])):
        off = int(q.split("|")[1])
        a = hits[q].get("ttna", (0, "", ""))
        b = hits[q].get("ttnb", (0, "", ""))
        t = "ttna" if a[0] >= b[0] else "ttnb"
        rows.append({"model": acc, "window_start_aa": off + 1, "best_zebrafish_paralogue": t,
                     "score": max(a[0], b[0]), "identity": "",
                     "score_ttna": a[0], "ttna_pos": f"{a[1]}-{a[2]}", "score_ttnb": b[0], "ttnb_pos": f"{b[1]}-{b[2]}"})
with open(os.path.join(D, "followups", "fish_titin_pair.tsv"), "w") as f:
    cols = ["model", "window_start_aa", "best_zebrafish_paralogue", "score_ttna", "ttna_pos", "score_ttnb", "ttnb_pos"]
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
for acc in ("KAK2904762.1", "TSR51413.1"):
    rs = [r for r in rows if r["model"] == acc]
    seq_order = []
    for r in rs:
        if not seq_order or seq_order[-1][0] != r["best_zebrafish_paralogue"]:
            seq_order.append([r["best_zebrafish_paralogue"], r["window_start_aa"], r["window_start_aa"], 1])
        else:
            seq_order[-1][2] = r["window_start_aa"]
            seq_order[-1][3] += 1
    print(acc, "windows with a hit:", len(rs))
    print("  runs of best paralogue along the model (paralogue, first window aa, last window aa, n windows):")
    for s in seq_order:
        print("   ", s)
