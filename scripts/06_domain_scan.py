"""Targeted, uniform domain scan of every unique giant sequence with pyhmmer.

A full Pfam scan of ~20 M residues is ~10-20 CPU-hours (estimate), outside the
Stage 0 budget; here we scan a fixed panel of ~70 Pfam HMMs that classify the
architectures expected among giant proteins (PKS, NRPS, titin/obscurin-like,
bacterial adhesin repeats, dynein, cytoskeletal linkers, SecA-like) and use
Pfam gathering thresholds. UniParc's pre-computed InterPro features are used
elsewhere for the full picture when they exist.

Input : data/fasta/unique_ge35992.fasta (written by 04_flags.py)
Output: data/hmm/panel.hmm, data/domain_panel_hits.tsv, data/domain_panel_summary.tsv
"""
import gzip
import io
import os
import sys
import time
from collections import Counter, defaultdict

import pyhmmer

sys.path.insert(0, os.path.dirname(__file__))
from http_util import get  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
os.makedirs(os.path.join(D, "hmm"), exist_ok=True)

PANEL = {
    # type I PKS
    "PF00109": "PKS", "PF02801": "PKS", "PF00698": "PKS", "PF08659": "PKS", "PF14765": "PKS",
    "PF08240": "PKS", "PF00107": "PKS", "PF00975": "PKS/NRPS", "PF00550": "carrier",
    "PF08242": "MT", "PF13649": "MT",
    # NRPS
    "PF00501": "NRPS", "PF13193": "NRPS", "PF00668": "NRPS",
    # titin / obscurin / muscle
    "PF07679": "Ig/FN3", "PF00041": "Ig/FN3", "PF13927": "Ig/FN3", "PF00047": "Ig/FN3",
    "PF13895": "Ig/FN3", "PF00069": "kinase", "PF09042": "titin", "PF00880": "nebulin",
    "PF00435": "spectrin", "PF00307": "CH", "PF00681": "plakin", "PF02187": "GAS2",
    "PF01576": "myosin", "PF00063": "myosin", "PF00612": "IQ",
    # dynein heavy chain
    "PF03028": "dynein", "PF12774": "dynein", "PF12775": "dynein", "PF12777": "dynein",
    "PF08393": "dynein", "PF12780": "dynein", "PF18198": "dynein", "PF18199": "dynein",
    # bacterial surface / adhesin repeats
    "PF05345": "adhesin", "PF00353": "adhesin", "PF13573": "adhesin", "PF02368": "adhesin",
    "PF07523": "adhesin", "PF00395": "adhesin", "PF00404": "adhesin", "PF00963": "adhesin",
    "PF13385": "adhesin", "PF00028": "adhesin", "PF07691": "adhesin", "PF13205": "adhesin",
    "PF01345": "adhesin", "PF17963": "adhesin", "PF07705": "adhesin", "PF05593": "adhesin",
    "PF03527": "adhesin", "PF13620": "adhesin",
    # SecA-like (Omnitrophota giants)
    "PF07517": "SecA", "PF01043": "SecA", "PF07516": "SecA",
    # misc eukaryotic
    "PF12796": "ankyrin", "PF00400": "WD40", "PF00092": "VWA", "PF00008": "EGF",
    "PF00057": "LDLa", "PF00090": "TSP1", "PF01391": "collagen", "PF00084": "sushi",
    "PF02026": "RyR", "PF00038": "filament", "PF00023": "ankyrin", "PF00515": "TPR",
    "PF13181": "TPR", "PF00560": "LRR", "PF13855": "LRR", "PF00004": "AAA",
}

panel_path = os.path.join(D, "hmm", "panel.hmm")
if not os.path.exists(panel_path):
    with open(panel_path, "wb") as out:
        for acc in PANEL:
            raw = get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}", {"annotation": "hmm"})
            try:
                raw = gzip.decompress(raw)
            except OSError:
                pass
            out.write(raw if raw.endswith(b"\n") else raw + b"\n")

with pyhmmer.plan7.HMMFile(panel_path) as hf:
    hmms = list(hf)
print("HMMs loaded:", len(hmms))
alphabet = pyhmmer.easel.Alphabet.amino()
# HMMER's pipeline refuses targets > 100,000 residues: cut every sequence into
# 40,000-aa windows overlapping by 1,500 aa and map hits back; a hit is kept only
# from the window where it starts inside the window's own (non-overlap) part.
WIN, OVL = 40000, 1500
with pyhmmer.easel.SequenceFile(os.path.join(D, "fasta", "unique_ge35992.fasta"), digital=True,
                                alphabet=alphabet) as sf:
    full = sf.read_block()
blocks = []
for sq in full:
    nm = sq.name.decode() if isinstance(sq.name, bytes) else sq.name
    text = sq.textize().sequence
    for off in range(0, len(text), WIN - OVL):
        piece = text[off:off + WIN]
        blocks.append(pyhmmer.easel.TextSequence(name=f"{nm}|{off}".encode(), sequence=piece).digitize(alphabet))
        if off + WIN >= len(text):
            break
seqs = pyhmmer.easel.DigitalSequenceBlock(alphabet, blocks)
def _s(x):
    return x.decode() if isinstance(x, bytes) else x


L = {(_s(s.name)): len(s) for s in full}
t0 = time.time()
cpu0 = time.process_time()
hits = []
for th in pyhmmer.hmmsearch(hmms, seqs, cpus=int(os.environ.get("HMM_CPUS", "8")), bit_cutoffs="gathering"):
    fam = _s(th.query.accession).split(".")[0]
    for h in th:
        if not h.included:
            continue
        nm, off = _s(h.name).rsplit("|", 1)
        off = int(off)
        for dmn in h.domains:
            if not dmn.included:
                continue
            if off > 0 and dmn.env_from <= OVL // 2:
                continue  # belongs to the previous window
            if dmn.env_from > WIN - OVL // 2 and off + WIN < L[nm]:
                continue  # belongs to the next window
            hits.append((nm, off + dmn.env_from, off + dmn.env_to, fam, PANEL.get(fam, "?"),
                         round(dmn.score, 1)))
wall = time.time() - t0
print(f"hmmsearch wall {wall:.1f}s, main-process cpu {time.process_time()-cpu0:.1f}s, residues {sum(L.values())}")

hits.sort()
with open(os.path.join(D, "domain_panel_hits.tsv"), "w") as f:
    f.write("md5\tstart\tend\tpfam\tclass\tscore\n")
    for h in hits:
        f.write("\t".join(map(str, h)) + "\n")

by = defaultdict(list)
for h in hits:
    by[h[0]].append(h)
with open(os.path.join(D, "domain_panel_summary.tsv"), "w") as f:
    f.write("md5\tlength\tpanel_hits\tpanel_cov_frac\tclass_counts\tpanel_class\tclass_order\n")
    for name, n in L.items():
        hs = by.get(name, [])
        cov = bytearray(n)
        for _, a, b, *_ in hs:
            cov[a - 1:b] = b"\x01" * (b - a + 1)
        cc = Counter(h[4] for h in hs)
        pks = cc["PKS"] >= 3 or (cc["PKS"] >= 2 and cc["carrier"] >= 1)
        nrps = cc["NRPS"] >= 2
        if pks and nrps:
            cls = "PKS+NRPS"
        elif pks:
            cls = "PKS"
        elif nrps:
            cls = "NRPS"
        elif cc["Ig/FN3"] >= 40:
            cls = "Ig/FN3-repeat (titin-like)" if cc["kinase"] or cc["titin"] else "Ig/FN3-repeat"
        elif cc["adhesin"] >= 10:
            cls = "bacterial-adhesin-repeat"
        elif cc["dynein"] >= 3:
            cls = "dynein-HC-like"
        elif cc:
            cls = "mixed/other:" + ",".join(k for k, _ in cc.most_common(3))
        else:
            cls = "no panel hit"
        # order of classes along the sequence, consecutive duplicates collapsed
        order = []
        for h in hs:
            if not order or order[-1][0] != h[4]:
                order.append([h[4], 1])
            else:
                order[-1][1] += 1
        f.write(f"{name}\t{n}\t{len(hs)}\t{sum(cov)/n:.3f}\t"
                f"{';'.join(f'{k}:{v}' for k, v in cc.most_common())}\t{cls}\t"
                f"{'>'.join(f'{k}*{v}' for k, v in order)[:300]}\n")
print("sequences with >=1 panel hit:", len(by), "of", len(L))
