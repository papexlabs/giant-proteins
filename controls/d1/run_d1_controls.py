"""Positive controls for the Stage 1 contradiction rules D1-D3 on real repetitive proteins.

SELECTION RULE (fixed on 13 Nov 2026 before any control was scored; queries run against the
UniProtKB REST API at run time, release in the output):
  G1  silks: reviewed (Swiss-Prot) entries, 'Evidence at protein level', length >= 300 aa,
      protein name containing fibroin, spidroin or silk, excluding 'silk gland factor'
      (transcription factors, not silk proteins)
  G1b complete spidroins: unreviewed entries, 'Evidence at protein level', not fragments,
      length >= 2,500 aa, protein name containing spidroin (Swiss-Prot holds almost no complete
      spidroin; these are real silk proteins with peptide evidence in UniProtKB)
  G2  mucins: reviewed, protein level, length >= 2,000 aa, protein name containing mucin
  G3  bacterial adhesins: reviewed, protein level, length >= 2,000 aa, Bacteria, protein name
      containing adhesin or keyword Cell adhesion (KW-0130)
  G4  long proteins: reviewed, protein level, length >= 5,000 aa (any organism)
An entry in several groups is scored once and reported in each group.

SCORING: exactly the Stage 1 rules of scripts/13_grade.py, with the same code paths:
  protein period and identity (scripts/period.py, lags 1-2,000, every 7th position), rep25
  (scripts/04_flags.py), domain coverage = max(panel coverage with data/hmm/panel.hmm at
  gathering thresholds, UniParc Pfam coverage from the entry's UniParc record), and, where an
  EMBL coding sequence is retrievable from ENA and translates to the entry's sequence (>= 95 %
  identical, same length +-1 %; at most 3 EMBL CDS tried per entry, mRNA-derived first),
  Fickett windows and open alternative frames as in
  scripts/11_locus_dna_checks.py.
  D1: (identity >= 0.6 or rep25 >= 0.8) and domain coverage < 5 %
  D2: >= 50 % of 300-bp windows with Fickett <= 0.74, and domain coverage < 5 %
  D3: >= 2 of 5 alternative frames with < 5 stops per 1,000 codons, and domain coverage < 5 %
Output: controls/d1/d1_controls.tsv, controls/d1/summary.json
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter

import pyhmmer

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from period import periods  # noqa: E402

OUT = os.path.dirname(os.path.abspath(__file__))
UA = "PapexResearchBot/1.0 (+https://papex.org)"
GROUPS = {
    "G1 silk (Swiss-Prot)": '(reviewed:true) AND (existence:1) AND (length:[300 TO *]) AND '
                            '((protein_name:fibroin) OR (protein_name:spidroin) OR (protein_name:silk)) '
                            'NOT (protein_name:"silk gland factor")',
    "G1b complete spidroins": '(reviewed:false) AND (existence:1) AND (fragment:false) AND (length:[2500 TO *]) AND '
                              '(protein_name:spidroin)',
    "G2 mucins": '(reviewed:true) AND (existence:1) AND (length:[2000 TO *]) AND (protein_name:mucin)',
    "G3 bacterial adhesins": '(reviewed:true) AND (existence:1) AND (length:[2000 TO *]) AND (taxonomy_id:2) AND '
                             '((protein_name:adhesin) OR (keyword:KW-0130))',
    "G4 long proteins": '(reviewed:true) AND (existence:1) AND (length:[5000 TO *])',
}


def get(url, tries=3):
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": UA}), timeout=60) as r:
                return r.read(), dict(r.headers)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(3 * (k + 1))


def search(q):
    url = ("https://rest.uniprot.org/uniprotkb/search?format=json&size=500&fields="
           "accession,length,protein_name,organism_name,protein_existence,fragment,xref_embl,sequence,uniparc_id"
           "&query=" + urllib.parse.quote(q))
    out, rel = [], None
    while url:
        body, hdr = get(url)
        rel = hdr.get("X-UniProt-Release") or hdr.get("x-uniprot-release") or rel
        out += json.loads(body)["results"]
        m = re.search(r'<([^>]+)>;\s*rel="next"', hdr.get("Link", "") or hdr.get("link", ""))
        url = m.group(1) if m else None
        time.sleep(0.5)
    return out, rel


def rep_frac(s, k=25):
    seen = Counter(s[i:i + k] for i in range(len(s) - k + 1))
    mask = bytearray(len(s))
    for i in range(len(s) - k + 1):
        if seen[s[i:i + k]] > 1:
            mask[i:i + k] = b"\x01" * k
    return sum(mask) / max(1, len(s))


# Fickett / frames: the definitions of scripts/11_locus_dna_checks.py, executed from its source
src = open(os.path.join(ROOT, "scripts", "11_locus_dna_checks.py")).read()
ns = {}
exec(src[src.index("COMP = str.maketrans"):src.index("def analyze(")], ns)
fickett, stops_per_k, COMP = ns["fickett"], ns["stops_per_k"], ns["COMP"]
CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(bases):
    for j, b in enumerate(bases):
        for k2, c in enumerate(bases):
            CODON[a + b + c] = aas[16 * i + 4 * j + k2]


def translate(x):
    return "".join(CODON.get(x[i:i + 3], "X") for i in range(0, len(x) - 2, 3))


entries, release = {}, None
member = {}
for g, q in GROUPS.items():
    res, rel = search(q)
    release = rel or release
    for e in res:
        acc = e["primaryAccession"]
        entries[acc] = e
        member.setdefault(acc, []).append(g)
print("entries", len(entries), "release", release, {g: sum(1 for a in member if g in member[a]) for g in GROUPS})

# panel coverage
alphabet = pyhmmer.easel.Alphabet.amino()
with pyhmmer.plan7.HMMFile(os.path.join(ROOT, "data", "hmm", "panel.hmm")) as hf:
    hmms = list(hf)
seqs = pyhmmer.easel.DigitalSequenceBlock(alphabet, [
    pyhmmer.easel.TextSequence(name=a.encode(), sequence=e["sequence"]["value"]).digitize(alphabet)
    for a, e in entries.items()])
cov = {a: bytearray(len(e["sequence"]["value"])) for a, e in entries.items()}
for th in pyhmmer.hmmsearch(hmms, seqs, cpus=4, bit_cutoffs="gathering"):
    for h in th:
        if not h.included:
            continue
        nm = h.name.decode() if isinstance(h.name, bytes) else h.name
        for dmn in h.domains:
            if dmn.included:
                cov[nm][dmn.env_from - 1:dmn.env_to] = b"\x01" * (dmn.env_to - dmn.env_from + 1)

rows = []
for n, (acc, e) in enumerate(sorted(entries.items(), key=lambda x: -x[1]["sequence"]["length"])):
    s = e["sequence"]["value"]
    L = len(s)
    top = periods(s, maxp=2000, top=1)[0]
    r25 = rep_frac(s)
    panel = sum(cov[acc]) / L
    upi = (e.get("extraAttributes") or {}).get("uniParcId", "")
    upf = None
    if upi:
        try:
            body, _ = get(f"https://rest.uniprot.org/uniparc/{upi}?fields=sequencefeatures&format=json")
            feats = json.loads(body).get("sequenceFeatures", [])
            pc = bytearray(L)
            for f in feats:
                if f.get("database") == "Pfam":
                    for loc in f.get("locations", []):
                        pc[loc["start"] - 1:min(L, loc["end"])] = b"\x01" * (min(L, loc["end"]) - loc["start"] + 1)
            upf = sum(pc) / L
        except Exception:
            upf = None
        time.sleep(0.4)
    dom = max(panel, upf or 0.0)
    # coding sequence from ENA (first EMBL protein_id whose translation matches)
    cds_info = {}
    cands = []
    for x in e.get("uniProtKBCrossReferences", []):
        if x.get("database") != "EMBL":
            continue
        props = {p["key"]: p["value"] for p in x.get("properties", [])}
        pid = props.get("ProteinId", "")
        if not pid or pid == "-" or props.get("Status", "-") not in ("-", ""):
            continue
        cands.append((props.get("MoleculeType") != "mRNA", pid))   # mRNA-derived CDS first
    for _, pid in sorted(cands)[:3]:                                  # at most 3 attempts per entry
        try:
            body, _ = get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{pid}")
        except Exception:
            continue
        time.sleep(0.4)
        cds = "".join(l.strip().upper() for l in body.decode().splitlines() if l and not l.startswith(">"))
        if len(cds) < 300:
            continue
        t = translate(cds).rstrip("*")
        ident = sum(1 for a, b in zip(t, s) if a == b) / max(len(t), L)
        if ident >= 0.95 and abs(len(t) - L) <= 0.01 * L + 1:
            fk = [fickett(cds[w:w + 300]) for w in range(0, len(cds) - 299, 300)]
            rc = cds.translate(COMP)[::-1]
            other = [stops_per_k(cds, 1), stops_per_k(cds, 2), stops_per_k(rc, 0), stops_per_k(rc, 1), stops_per_k(rc, 2)]
            cds_info = {"cds_protein_id": pid, "cds_translation_identity": round(ident, 3),
                        "fickett_noncoding_frac": round(sum(1 for v in fk if v <= 0.74) / len(fk), 3),
                        "fickett_median": round(sorted(fk)[len(fk) // 2], 3),
                        "n_other_frames_open": sum(1 for v in other if v < 5),
                        "min_stops_per_kcodon_other_frames": round(min(other), 2)}
            break
    d1 = (top[0] >= 0.6 or r25 >= 0.8) and dom < 0.05
    d2 = bool(cds_info) and cds_info["fickett_noncoding_frac"] >= 0.5 and dom < 0.05
    d3 = bool(cds_info) and cds_info["n_other_frames_open"] >= 2 and dom < 0.05
    rows.append({
        "accession": acc, "groups": "; ".join(member[acc]), "length": L,
        "protein_name": ((e.get("proteinDescription") or {}).get("recommendedName") or
                         ((e.get("proteinDescription") or {}).get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", ""),
        "organism": (e.get("organism") or {}).get("scientificName", ""),
        "existence": e.get("proteinExistence", ""), "fragment": bool((e.get("sequence") or {}).get("fragment") or
                                                                      "fragment" in json.dumps(e.get("proteinDescription", {})).lower()),
        "period_aa": top[1], "period_identity": round(top[0], 3), "rep25": round(r25, 3),
        "panel_cov": round(panel, 3), "uniparc_pfam_cov": "" if upf is None else round(upf, 3),
        "domain_cov": round(dom, 3), "D1": d1, "D2": d2, "D3": d3,
        "d1_trigger_without_domain_guard": top[0] >= 0.6 or r25 >= 0.8, **cds_info})
    if n % 10 == 0:
        print(n, acc, L, rows[-1]["period_identity"], rows[-1]["rep25"], rows[-1]["domain_cov"], d1, flush=True)

cols = []
for r in rows:
    for k in r:
        if k not in cols:
            cols.append(k)
with open(os.path.join(OUT, "d1_controls.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
summ = {"uniprot_release": release, "n_entries": len(rows), "queries": GROUPS, "by_group": {}}
for g in GROUPS:
    rr = [r for r in rows if g in r["groups"]]
    summ["by_group"][g] = {"n": len(rr), "D1": sum(r["D1"] for r in rr), "D2": sum(r["D2"] for r in rr),
                           "D3": sum(r["D3"] for r in rr), "any_D": sum(r["D1"] or r["D2"] or r["D3"] for r in rr),
                           "with_cds": sum(1 for r in rr if r.get("cds_protein_id")),
                           "repeat_trigger_before_domain_guard": sum(r["d1_trigger_without_domain_guard"] for r in rr)}
summ["all"] = {"n": len(rows), "D1": sum(r["D1"] for r in rows), "any_D": sum(r["D1"] or r["D2"] or r["D3"] for r in rows),
               "flagged": [(r["accession"], r["protein_name"][:50], r["length"], r["period_aa"], r["period_identity"],
                            r["rep25"], r["domain_cov"], r["D1"], r["D2"], r["D3"]) for r in rows if r["D1"] or r["D2"] or r["D3"]]}
json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
print(json.dumps(summ["by_group"], indent=1))
print("flagged:", summ["all"]["flagged"])
