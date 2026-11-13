"""Follow-up (c): why do five titin models have low junction support (B -> C)?

For each model and run (data from the RNA-seq stage results):
  - library metadata from ENA (strategy, selection, layout, read length, instrument)
  - junction fraction: reads on annotated junctions / unique reads (human titin control: 0.172)
  - 5'->3' profile in five equal bins of the exon chain (transcription order): mean exon depth
    and share of introns with >= 3 exact reads; a 3'-tag or degraded library piles up at the
    3' end, a non-expressing tissue is low everywhere
  - "covered but unspliced" introns: both flanking exons at depth >= 10 yet 0 exact reads;
    among those, introns whose donor or acceptor IS used with another partner (alternative
    junction reads): that points at a wrong exon boundary in the model rather than at the run
The human titin control is processed the same way for comparison.
Output: data/followups/titin_downgrades.tsv and printed summary.
"""
import csv
import glob
import json
import os
import urllib.parse
import urllib.request

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
UA = "PapexResearchBot/1.0 (+https://papex.org)"
MODELS = {"XP_040214433.1": "Rana temporaria", "XP_030066708.1": "Microcaecilia unicolor",
          "XP_071978768.1": "Engystomops pustulosus", "XP_087738184.1": "Pantherophis guttatus",
          "KAG7502825.1": "Solea senegalensis", "NP_001254479.2": "Homo sapiens (control)"}
loci = {}
for p in glob.glob(os.path.join(ROOT, "cloud", "bundle", "loci", "*.json")):
    for L in json.load(open(p)):
        loci[L["rep_accession"]] = L


def ena_meta(run):
    q = urllib.parse.urlencode({"result": "read_run", "query": f'run_accession="{run}"', "format": "tsv",
                                "fields": "library_strategy,library_selection,library_source,library_layout,"
                                          "instrument_model,read_count,base_count,sample_title"})
    t = urllib.request.urlopen(urllib.request.Request("https://www.ebi.ac.uk/ena/portal/api/search?" + q,
                                                      headers={"User-Agent": UA}), timeout=60).read().decode()
    lines = t.splitlines()
    if len(lines) < 2:
        return {}
    r = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    try:
        mates = 2 if r.get("library_layout") == "PAIRED" else 1
        r["read_len"] = round(int(r["base_count"]) / int(r["read_count"]) / mates)
    except (ValueError, ZeroDivisionError, KeyError):
        r["read_len"] = ""
    return r


rows = []
for acc, org in MODELS.items():
    L0 = loci.get(acc)
    for p in sorted(glob.glob(os.path.join(ROOT, "cloud_out", "results", "*", "*.json"))):
        if p.endswith(".skip.json"):
            continue
        d = json.load(open(p))
        L = next((x for x in d.get("loci", []) if x["rep_accession"] == acc), None)
        if not L:
            continue
        strand = (L0 or {}).get("strand", "+")
        ex = sorted(L["exons"], key=lambda e: e[0], reverse=(strand == "-"))     # transcription order
        intr = sorted(L["introns"], key=lambda i: i[0], reverse=(strand == "-"))
        n_e, n_i = len(ex), len(intr)
        bins = []
        for b in range(5):
            es = ex[b * n_e // 5:(b + 1) * n_e // 5]
            is_ = intr[b * n_i // 5:(b + 1) * n_i // 5]
            bins.append((round(sum(e[3] for e in es) / max(1, len(es)), 1),
                         round(sum(1 for i in is_ if i[2] >= 3) / max(1, len(is_)), 2)))
        depth_at = {}
        for e in L["exons"]:
            depth_at[e[1]] = e[3]          # exon ending at a donor (genomic left exon)
            depth_at[("s", e[0])] = e[3]   # exon starting at an acceptor
        covered_unspliced = alt_used = 0
        for a, b, u, m, ad, aa in L["introns"]:
            left = depth_at.get(a - 1, 0)
            right = depth_at.get(("s", b + 1), 0)
            if left >= 10 and right >= 10 and u + m == 0:
                covered_unspliced += 1
                if ad + aa > 0:
                    alt_used += 1
        ann = sum(i[2] + i[3] for i in L["introns"])
        meta = ena_meta(d["run"])
        rows.append({"model": acc, "organism": org, "run": d["run"],
                     "sample": (d.get("run_meta", {}).get("sample") or "")[:40],
                     "library": f"{meta.get('library_strategy','')}/{meta.get('library_selection','')}/"
                                f"{meta.get('library_layout','')}",
                     "read_len": meta.get("read_len", ""), "instrument": meta.get("instrument_model", ""),
                     "reads_unique": L["reads_unique"], "junction_fraction": round(ann / max(1, L["reads_unique"]), 3),
                     "introns_ge3": L["introns_with_ge3_exact_reads"], "n_introns": n_i,
                     "profile_5to3_depth": [b[0] for b in bins], "profile_5to3_supported": [b[1] for b in bins],
                     "covered_but_unspliced_introns": covered_unspliced,
                     "of_which_partner_used_elsewhere": alt_used})
cols = list(rows[0].keys())
with open(os.path.join(D, "followups", "titin_downgrades.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")
for r in rows:
    print(f"{r['organism'][:22]:22s} {r['run']:12s} {r['sample'][:24]:24s} {r['library'][:28]:28s} len {r['read_len']!s:>4} "
          f"uniq {r['reads_unique']:>7} jfrac {r['junction_fraction']:<5} int {r['introns_ge3']:>3}/{r['n_introns']} "
          f"5'->3' depth {r['profile_5to3_depth']} supp {r['profile_5to3_supported']} "
          f"cov-unspliced {r['covered_but_unspliced_introns']} (alt {r['of_which_partner_used_elsewhere']})")
