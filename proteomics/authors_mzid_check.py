"""Where do the authors' own identifications on the giant fall? (one of their six Scaffold mzIdentML exports)

Usage: python3 authors_mzid_check.py raw/mzid/<file>.mzid.gz
Counts PeptideEvidence entries on DBSeq_Omnitrophica_OP3_32 (the giant in the authors' 2017 database),
the distinct peptides behind the identification items that reference them, their start positions per
fifth of the length, and how many of those items carry passThreshold="true".
Writes results/authors_mzid_check.json. Streaming XML parse (iterparse).
"""
import gzip
import json
import os
import sys
import xml.etree.ElementTree as ET
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
NS = "{http://psidev.info/psi/pi/mzIdentML/1.1}"
GI = "DBSeq_Omnitrophica_OP3_32"
L, ev, items, passed, peps_hit = None, {}, 0, 0, set()
starts_hit = set()
starts_pass, peps_pass = set(), set()
for _, el in ET.iterparse(gzip.open(sys.argv[1]), events=("end",)):
    t = el.tag
    if t == NS + "DBSequence":
        if el.get("id") == GI:
            L = int(el.get("length"))
        el.clear()
    elif t == NS + "PeptideEvidence":
        if el.get("dBSequence_ref") == GI:
            ev[el.get("id")] = (int(el.get("start")), el.get("peptide_ref"), el.get("isDecoy"))
        el.clear()
    elif t == NS + "SpectrumIdentificationItem":
        refs = [r.get("peptideEvidence_ref") for r in el.findall(NS + "PeptideEvidenceRef")]
        on = [r for r in refs if r in ev]
        if on:
            items += 1
            if el.get("passThreshold") == "true":
                passed += 1
                for r in on:
                    starts_pass.add(ev[r][0])
                    peps_pass.add(ev[r][1])
            for r in on:
                starts_hit.add(ev[r][0])
                peps_hit.add(ev[r][1])
        el.clear()
    elif t == NS + "SpectrumIdentificationResult":
        el.clear()
starts = sorted(starts_hit)
fifth = Counter(min(4, (s - 1) * 5 // L) for s in starts)
out = {"file": os.path.basename(sys.argv[1]), "giant_length_in_authors_db": L,
       "peptide_evidence_on_giant": len(ev), "decoy_evidence": sum(1 for v in ev.values() if v[2] == "true"),
       "identification_items_on_giant": items, "items_passThreshold_true": passed,
       "distinct_peptides": len(peps_hit), "distinct_start_positions": len(starts),
       "starts_per_fifth": [fifth[i] for i in range(5)], "first_start": starts[0] if starts else None,
       "last_start": starts[-1] if starts else None, "starts_after_6000": [s for s in starts if s > 6000],
       "passThreshold_only": {"distinct_peptides": len(peps_pass), "distinct_start_positions": len(starts_pass),
                              "starts_per_fifth": [Counter(min(4, (s - 1) * 5 // L) for s in starts_pass)[i] for i in range(5)],
                              "starts_after_6000": sorted(s for s in starts_pass if s > 6000)}}
json.dump(out, open(os.path.join(HERE, "results", "authors_mzid_check.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
