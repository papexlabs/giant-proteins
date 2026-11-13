#!/usr/bin/env python3
"""How many introns of each repeat-flag model can get unique junction support at all?

For every model in cloud/bundle_rep_rna/loci/*.json, the spliced junction sequence of each
annotated intron (the last 50 bp of the upstream exon + the first 50 bp of the downstream exon,
in transcript orientation; shorter exons contribute what they have) is compared with the
junction sequences of all other introns of the same model and with the model's own exon
sequence (both strands). An intron whose junction 100-mer occurs elsewhere in the model cannot
be supported by a uniquely placed read (MAPQ >= 10), because a read across it fits at least
two places in the gene. This is an upper bound: copies elsewhere in the genome are not checked.

Also counted: introns whose two flanking exons both contain a 50-mer unique within the model
("flanks mappable"). Those can become testable (unique flank depth >= 10) even when their
junction is not unique: the configuration in which a correct model could be scored as missing
support.

Output: controls/rna_rep/junction_mappability.tsv and .json (written BEFORE any RNA-seq result
of these models exists; part of the pre-registered rule, see cloud/bundle_rep_rna/RULES.txt).
Usage: python3 controls/rna_rep/junction_mappability.py
"""
import glob
import gzip
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
D = os.path.join(ROOT, "data")
OUT = os.path.dirname(os.path.abspath(__file__))
IDX = json.load(open(os.path.join(D, "loci", "dna_index.json")))
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_dna = {}


def region(md5):
    for lid, v in IDX.items():
        if md5 in v["md5s"] and v["how"] not in ("fail", "fail-too-large"):
            if lid not in _dna:
                _dna[lid] = "".join(l.strip() for l in gzip.open(os.path.join(D, "loci", "dna", lid + ".fa.gz"), "rt")
                                    if not l.startswith(">")).upper()
            return _dna[lid], v["start"]
    return "", 0


rows = []
for f in sorted(glob.glob(os.path.join(ROOT, "cloud", "bundle_rep_rna", "loci", "*.json"))):
    sp = os.path.basename(f)[:-5]
    for L in json.load(open(f)):
        seq, start = region(L["md5"])
        if not seq:
            rows.append({"species_key": sp, "rep_accession": L["rep_accession"], "n_introns": "", "note": "no locus DNA"})
            continue
        segs = sorted(L["segments"])                      # genomic order, + strand coordinates
        ex = [seq[a - start:b - start + 1] for a, b in segs]
        if L["strand"] == "-":
            ex = [e.translate(COMP)[::-1] for e in ex][::-1]   # transcript orientation
        juncs = [ex[i][-50:] + ex[i + 1][:50] for i in range(len(ex) - 1)]
        cds = "".join(ex)
        rcds = cds.translate(COMP)[::-1]

        def occ(j, s, stop=2):
            n, i = 0, s.find(j)
            while i >= 0 and n < stop:
                n += 1
                i = s.find(j, i + 1)
            return n
        # the transcript (cds) contains every junction once at its own place; unique = no second
        # (also overlapping) occurrence on the same strand and none on the other strand
        uj = [occ(j, cds) == 1 and occ(j, rcds, 1) == 0 for j in juncs]
        uniq_j = sum(uj)
        k50 = {}
        for i in range(len(cds) - 49):
            k = cds[i:i + 50]
            k50[k] = k50.get(k, 0) + 1
        for i in range(len(rcds) - 49):
            k = rcds[i:i + 50]
            if k in k50:
                k50[k] += 1

        def has_unique(e):
            return any(k50.get(e[i:i + 50], 0) == 1 for i in range(max(0, len(e) - 49)))
        flank_ok = [has_unique(e) for e in ex]
        flanks = sum(1 for i in range(len(ex) - 1) if flank_ok[i] and flank_ok[i + 1])
        flanks_not_junc = sum(1 for i in range(len(juncs)) if flank_ok[i] and flank_ok[i + 1] and not uj[i])
        n = len(juncs)
        # per intron in GENOMIC order (the order of rnaseq_locus.py's intron list): [length, unique junction]
        g_uj = uj[::-1] if L["strand"] == "-" else uj
        per = [[segs[i + 1][0] - segs[i][1] - 1, int(g_uj[i])] for i in range(n)]
        rows.append({"species_key": sp, "rep_accession": L["rep_accession"], "md5": L["md5"], "n_introns": n,
                     "per_intron_genomic_order": per,
                     "introns_unique_junction": uniq_j, "share_unique_junction": round(uniq_j / n, 3) if n else "",
                     "introns_flanks_mappable": flanks, "flanks_mappable_junction_not_unique": flanks_not_junc,
                     "note": ""})
cols = ["species_key", "rep_accession", "n_introns", "introns_unique_junction", "share_unique_junction",
        "introns_flanks_mappable", "flanks_mappable_junction_not_unique", "note"]
with open(os.path.join(OUT, "junction_mappability.tsv"), "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in rows:
        fo.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
json.dump(rows, open(os.path.join(OUT, "junction_mappability.json"), "w"), indent=1)
for r in rows:
    print(r["rep_accession"], r["species_key"][:30], r.get("n_introns"), r.get("introns_unique_junction"),
          r.get("share_unique_junction"), r.get("introns_flanks_mappable"), r.get("flanks_mappable_junction_not_unique"), r["note"])
