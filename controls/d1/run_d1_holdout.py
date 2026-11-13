"""Independent hold-out test of the repeat rule (D1) and of its revision (D1 counts as a
contradiction only together with D2 or D3).

SELECTION RULE (fixed on 13 Nov 2026 after the design controls in run_d1_controls.py had been
scored, and before this set was scored): every UniProtKB/Swiss-Prot entry with 'Evidence at
protein level' and length 1,000-4,999 aa that is not already in the design set
(controls/d1/d1_controls.tsv). These are real proteins of every kind; the repetitive ones
among them are the test.
SCORING: the protein-level D1 trigger exactly as scripts/13_grade.py (period identity >= 0.6
at the best lag 1-2,000, every 7th position; or rep25 >= 0.8), the domain guard (max of the
74-HMM panel coverage and UniParc Pfam coverage < 5 %), and for the triggered entries the
coding sequence from ENA (at most 3 EMBL CDS, mRNA first, translation >= 95 % identical) for
D2 (>= 50 % non-coding Fickett windows) and D3 (>= 2 of 5 alternative frames open).
Output: controls/d1/holdout.tsv (triggered entries only) and holdout_summary.json
"""
import csv
import json
import os
import sys
import time
import urllib.parse

import numpy as np
import pyhmmer

HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = [sys.argv[0]]
src = open(os.path.join(HERE, "run_d1_controls.py")).read()
ns = {"__file__": os.path.join(HERE, "run_d1_controls.py")}
exec(src[:src.index("entries, release = {}, None")], ns)          # helpers only (get, rep_frac, fickett, ...)
get, rep_frac, fickett, stops_per_k, COMP, translate, UA, ROOT = (ns[k] for k in (
    "get", "rep_frac", "fickett", "stops_per_k", "COMP", "translate", "UA", "ROOT"))

design = {r["accession"] for r in csv.DictReader(open(os.path.join(HERE, "d1_controls.tsv")), delimiter="\t")}
q = "(reviewed:true) AND (existence:1) AND (length:[1000 TO 4999])"
body, hdr = get("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=" + urllib.parse.quote(q))
release = hdr.get("X-UniProt-Release") or hdr.get("x-uniprot-release")
seqs, name = {}, None
for line in body.decode().splitlines():
    if line.startswith(">"):
        name = line.split("|")[1]
        seqs[name] = []
    elif name:
        seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items() if k not in design}
print("hold-out entries", len(seqs), "release", release, flush=True)


def period_np(s, maxp=2000):
    a = np.frombuffer(s.encode(), dtype=np.uint8)
    n = len(a)
    best = (0.0, 0)
    for p in range(1, min(maxp, n // 2)):
        idx = np.arange(0, n - p, 7)
        m = float((a[idx] == a[idx + p]).mean())
        if m > best[0]:
            best = (m, p)
    return best


trig = {}
t0 = time.time()
for i, (acc, s) in enumerate(seqs.items()):
    pid, per = period_np(s)
    r25 = rep_frac(s)
    if pid >= 0.6 or r25 >= 0.8:
        trig[acc] = {"period_aa": per, "period_identity": round(pid, 3), "rep25": round(r25, 3)}
    if i % 1000 == 0:
        print(i, len(trig), f"{time.time() - t0:.0f}s", flush=True)
print("repeat trigger (before the domain guard):", len(trig), flush=True)

alphabet = pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile(os.path.join(ROOT, "data", "hmm", "panel.hmm")) as hf:
    hmms = list(hf)
blk = pyhmmer.easel.DigitalSequenceBlock(alphabet, [
    pyhmmer.easel.TextSequence(name=a.encode(), sequence=seqs[a]).digitize(alphabet) for a in trig])
cov = {a: bytearray(len(seqs[a])) for a in trig}
if trig:
    for th in pyhmmer.hmmsearch(hmms, blk, cpus=4, bit_cutoffs="gathering"):
        for h in th:
            if h.included:
                nm = h.name.decode() if isinstance(h.name, bytes) else h.name
                for d in h.domains:
                    if d.included:
                        cov[nm][d.env_from - 1:d.env_to] = b"\x01" * (d.env_to - d.env_from + 1)
rows = []
for acc, t in trig.items():
    s = seqs[acc]
    L = len(s)
    panel = sum(cov[acc]) / L
    body, _ = get(f"https://rest.uniprot.org/uniprotkb/{acc}?format=json&fields=protein_name,organism_name,xref_embl,uniparc_id")
    e = json.loads(body)
    time.sleep(0.3)
    upf = 0.0
    upi = (e.get("extraAttributes") or {}).get("uniParcId", "")
    if upi:
        try:
            b2, _ = get(f"https://rest.uniprot.org/uniparc/{upi}?fields=sequencefeatures&format=json")
            pc = bytearray(L)
            for f in json.loads(b2).get("sequenceFeatures", []):
                if f.get("database") == "Pfam":
                    for loc in f.get("locations", []):
                        pc[loc["start"] - 1:min(L, loc["end"])] = b"\x01" * (min(L, loc["end"]) - loc["start"] + 1)
            upf = sum(pc) / L
        except Exception:
            pass
        time.sleep(0.3)
    dom = max(panel, upf)
    row = {"accession": acc, "length": L,
           "protein_name": ((e.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName", {}).get("value", ""),
           "organism": (e.get("organism") or {}).get("scientificName", ""), **t, "domain_cov": round(dom, 3),
           "D1": dom < 0.05}
    if row["D1"]:
        cands = []
        for x in e.get("uniProtKBCrossReferences", []):
            if x.get("database") != "EMBL":
                continue
            props = {p["key"]: p["value"] for p in x.get("properties", [])}
            pid = props.get("ProteinId", "")
            if pid and pid != "-" and props.get("Status", "-") in ("-", ""):
                cands.append((props.get("MoleculeType") != "mRNA", pid))
        for _, pid in sorted(cands)[:3]:
            try:
                b3, _ = get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{pid}")
            except Exception:
                continue
            time.sleep(0.3)
            cds = "".join(l.strip().upper() for l in b3.decode().splitlines() if l and not l.startswith(">"))
            if len(cds) < 300:
                continue
            tr = translate(cds).rstrip("*")
            ident = sum(1 for a_, b_ in zip(tr, s) if a_ == b_) / max(len(tr), L)
            if ident >= 0.95 and abs(len(tr) - L) <= 0.01 * L + 1:
                fk = [fickett(cds[w:w + 300]) for w in range(0, len(cds) - 299, 300)]
                rc = cds.translate(COMP)[::-1]
                other = [stops_per_k(cds, 1), stops_per_k(cds, 2), stops_per_k(rc, 0), stops_per_k(rc, 1), stops_per_k(rc, 2)]
                row.update({"cds_protein_id": pid, "fickett_noncoding_frac": round(sum(1 for v in fk if v <= 0.74) / len(fk), 3),
                            "fickett_median": round(sorted(fk)[len(fk) // 2], 3),
                            "n_other_frames_open": sum(1 for v in other if v < 5)})
                break
        row["D2"] = bool(row.get("cds_protein_id")) and row["fickett_noncoding_frac"] >= 0.5
        row["D3"] = bool(row.get("cds_protein_id")) and row["n_other_frames_open"] >= 2
        row["D1_revised"] = row["D2"] or row["D3"]
    rows.append(row)
cols = []
for r in rows:
    for k in r:
        if k not in cols:
            cols.append(k)
with open(os.path.join(HERE, "holdout.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
d1 = [r for r in rows if r["D1"]]
summ = {"uniprot_release": release, "n_holdout": len(seqs), "repeat_trigger": len(trig), "D1_rule_v1": len(d1),
        "D1_with_cds": sum(1 for r in d1 if r.get("cds_protein_id")),
        "D1_revised_D1_and_D2_or_D3": sum(1 for r in d1 if r.get("D1_revised")),
        "D1_flagged": [(r["accession"], r["protein_name"][:50], r["organism"][:30], r["length"], r["period_aa"],
                        r["period_identity"], r["rep25"], r.get("fickett_median"), r.get("n_other_frames_open"),
                        r.get("D1_revised")) for r in d1]}
json.dump(summ, open(os.path.join(HERE, "holdout_summary.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != "D1_flagged"}))
for x in summ["D1_flagged"]:
    print(x)
