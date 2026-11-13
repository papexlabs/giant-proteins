"""Window homology: cut a query giant into fixed windows and search each window
with phmmer against a chosen set of giant sequences (by organism substring).
Reports, per target, how many query windows hit it (E <= 1e-10) and the
target positions, which shows whether two giants are homologous along their
whole length or only in parts (fusion / partial similarity).

Usage: window_homology.py <query rep_accession> <organism substring> [window=2000]
"""
import csv
import os
import sys

import pyhmmer

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
q_acc, org_sub = sys.argv[1], sys.argv[2]
W = int(sys.argv[3]) if len(sys.argv) > 3 else 2000

meta = {r["md5"]: r for r in csv.DictReader(open(os.path.join(D, "stage0_table.tsv")), delimiter="\t")}
alpha = pyhmmer.easel.Alphabet.amino()
seqs = {}
with pyhmmer.easel.SequenceFile(os.path.join(D, "fasta", "unique_ge35992.fasta"), digital=False) as sf:
    for s in sf:
        name = s.name.decode() if isinstance(s.name, bytes) else s.name
        seqs[name] = s.sequence
q_md5 = next(m for m, r in meta.items() if r["rep_accession"] == q_acc)
targets = [m for m, r in meta.items() if org_sub.lower() in r["organism"].lower() and m != q_md5
           and r["live"] == "yes"]
tblock = pyhmmer.easel.DigitalSequenceBlock(alpha, [
    pyhmmer.easel.TextSequence(name=m.encode(), sequence=seqs[m]).digitize(alpha) for m in targets])
qs = seqs[q_md5]
queries = [pyhmmer.easel.TextSequence(name=f"w{off}".encode(), sequence=qs[off:off + W]).digitize(alpha)
           for off in range(0, len(qs), W) if len(qs[off:off + W]) >= 200]
hitmap = {}
for th in pyhmmer.hmmer.phmmer(queries, tblock, cpus=8, E=1e-10):
    qn = th.query.name.decode() if isinstance(th.query.name, bytes) else th.query.name
    for h in th:
        if h.evalue > 1e-10:
            continue
        tn = h.name.decode() if isinstance(h.name, bytes) else h.name
        best = max(h.domains, key=lambda d: d.score)
        hitmap.setdefault(tn, []).append((int(qn[1:]), best.env_from, best.env_to, round(h.evalue, 30)))
print(f"query {q_acc} ({len(qs)} aa) in {len(queries)} windows of {W} aa; targets: {len(targets)}")
for tn, hs in sorted(hitmap.items(), key=lambda kv: -len(kv[1])):
    r = meta[tn]
    hs.sort()
    print(f"  {r['rep_accession']:16s} {r['length']:>6s} aa  {r['organism'][:32]:32s} windows hit {len(hs):>3d}/{len(queries)}"
          f"  query windows {hs[0][0]}..{hs[-1][0]}  target span {min(h[1] for h in hs)}..{max(h[2] for h in hs)}")
missing = [m for m in targets if m not in hitmap]
print(f"  targets with no window hit: {len(missing)}")
