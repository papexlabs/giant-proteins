"""Job 1, step 1: resolve the genomic locus of every non-excluded giant sequence
and download the locus DNA (span + 5 kb flanks), largest sequence first.

Locus sources, in order of preference:
  INSDC protein with a genomic /coded_by          -> exon chain from coded_by
  RefSeq XP_ (coded_by points to an XM_ mRNA)     -> NCBI Datasets gene product report
                                                     (transcript exons + CDS range)
  RefSeq WP_ (non-redundant, no locus)            -> NCBI identical-protein report (IPG)
  UniParc-only EMBL/EMBLWGS proteins              -> ENA EMBL record CDS location
DNA: NCBI efetch (seq_start/seq_stop); for contigs NCBI does not serve (e.g. ENA
TPA MAG sets) the whole ENA WGS set FASTA is downloaded if it is < 60 MB.
Resumable: specs in data/loci/specs.jsonl, DNA in data/loci/dna/<locus>.fa.gz.
No API keys or contact parameters; NCBI <= 2.5 req/s (http_util).
Usage: 10_locus_fetch.py [max_sequences]
"""
import csv
import gzip
import json
import os
import re
import sys
import urllib.error

sys.path.insert(0, os.path.dirname(__file__))
from http_util import get  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
LD = os.path.join(D, "loci")
os.makedirs(os.path.join(LD, "dna"), exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
FLANK = 5000
MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 9

stage0 = [r for r in csv.DictReader(open(os.path.join(D, "stage0_table.tsv")), delimiter="\t")
          if r["triage"] != "T-excl"]
stage0.sort(key=lambda r: -int(r["length"]))
rec = {r["accession"]: r for r in csv.DictReader(open(os.path.join(D, "ncbi_records_parsed.tsv")), delimiter="\t")}
xrefs = json.load(open(os.path.join(D, "raw", "uniparc_xrefs.json")))
# full /coded_by strings (the parsed TSV truncates them)
CODED = {}
for gp in gzip.open(os.path.join(D, "raw", "ncbi_gp.txt.gz"), "rt").read().split("\n//\n"):
    mv = re.search(r"^VERSION\s+(\S+)", gp, re.M)
    mc = re.search(r'/coded_by="(.*?)"', gp, re.S)
    if mv and mc:
        CODED[mv.group(1)] = re.sub(r"\s+", "", mc.group(1))

SPECS = os.path.join(LD, "specs.jsonl")
done = {}
if os.path.exists(SPECS):
    for line in open(SPECS):
        s = json.loads(line)
        done[s["md5"]] = s


def parse_location(loc):
    """INSDC location -> (contig, strand, [(start, end), ...] in transcription order)."""
    loc = re.sub(r"\s+", "", loc)
    whole_minus = loc.startswith("complement(") and loc.endswith(")") and \
        not re.match(r"^complement\([^()]*\),", loc)
    parts = re.findall(r"(complement\()?([A-Za-z0-9_]+\.\d+):<?(\d+)\.\.>?(\d+)", loc)
    if not parts:
        return None
    contigs = {p[1] for p in parts}
    segs = [(int(p[2]), int(p[3])) for p in parts]
    if whole_minus:
        return sorted(contigs)[0], "-", list(reversed(segs)), len(contigs)
    strand = "-" if parts[0][0] else "+"
    return sorted(contigs)[0], strand, segs, len(contigs)


gene_cache = {}


def from_gene(gene_id, xp):
    if gene_id not in gene_cache:
        gene_cache[gene_id] = json.loads(get(f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}/product_report"))
    for rep in gene_cache[gene_id].get("reports", []):
        for t in rep["product"].get("transcripts", []):
            if t.get("protein", {}).get("accession_version") != xp:
                continue
            gl = t["genomic_locations"][0]
            ex = sorted(gl["exons"], key=lambda e: int(e["order"]))
            cds = t["cds"]["range"][0]
            cb, ce = int(cds["begin"]), int(cds["end"])
            strand = "+" if gl["genomic_range"]["orientation"] == "plus" else "-"
            segs, pos = [], 0
            for e in ex:
                b, e_ = int(e["begin"]), int(e["end"])
                L = e_ - b + 1
                lo, hi = max(cb, pos + 1), min(ce, pos + L)        # transcript coords inside this exon
                if lo <= hi:
                    if strand == "+":
                        segs.append((b + lo - pos - 1, b + hi - pos - 1))
                    else:
                        segs.append((e_ - (lo - pos - 1), e_ - (hi - pos - 1)))
                pos += L
            # for '-' strand keep (start<=end) tuples
            segs = [(min(a, b), max(a, b)) for a, b in segs]
            return gl["genomic_accession_version"], strand, segs, 1, t["length"] - pos
    return None


def from_ipg(wp):
    t = get(EU + "efetch.fcgi", {"db": "protein", "id": wp, "rettype": "ipg", "retmode": "text"}).decode()
    for line in t.splitlines()[1:]:
        f = line.split("\t")
        if len(f) > 6 and f[2] and f[3].isdigit():
            return f[2], f[5], [(int(f[3]), int(f[4]))], 1
    return None


def from_ena(pid):
    t = get(f"https://www.ebi.ac.uk/ena/browser/api/embl/{pid}").decode()
    m = re.search(r"^FT   CDS\s+(.*?)\nFT {19}/", t, re.S | re.M)
    if not m:
        return None
    loc = re.sub(r"\nFT\s+", "", m.group(1))
    return parse_location(loc)


def contig_lengths(accs):
    out = {}
    accs = [a for a in accs if a]
    uids = {}
    for a in accs:                       # accession -> uid (esummary rejects some accessions)
        try:
            r = json.loads(get(EU + "esearch.fcgi", {"db": "nuccore", "term": f"{a}[accn]", "retmode": "json"}))
            if r["esearchresult"]["idlist"]:
                uids[r["esearchresult"]["idlist"][0]] = a
        except Exception:  # noqa: BLE001
            pass
    accs = list(uids)
    for i in range(0, len(accs), 100):
        try:
            r = json.loads(get(EU + "esummary.fcgi", data={"db": "nuccore", "id": ",".join(accs[i:i + 100]),
                                                           "retmode": "json"}))
            for uid in r["result"].get("uids", []):
                d = r["result"][uid]
                out[d["accessionversion"]] = int(d["slen"])
        except Exception:  # noqa: BLE001
            pass
    return out


def fetch_dna(contig, start, end, path):
    try:
        t = get(EU + "efetch.fcgi", {"db": "nuccore", "id": contig, "rettype": "fasta", "retmode": "text",
                                     "seq_start": start, "seq_stop": end}, timeout=120).decode()
        if t.startswith(">") and len(t) > 100:
            with gzip.open(path, "wt") as f:
                f.write(t)
            return "ncbi"
    except urllib.error.HTTPError:
        pass
    # ENA WGS set fallback (small sets only, e.g. MAGs)
    m = re.match(r"^([A-Z]{4,6})(\d{2})\d+\.\d+$", contig)
    if not m:
        return "fail"
    setname = m.group(1) + m.group(2)
    url = f"https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/{m.group(1)[:3].lower()}/{setname}.fasta.gz"
    setpath = os.path.join(LD, "dna", f"ENA_{setname}.fasta.gz")
    if not os.path.exists(setpath):
        try:
            raw = get(url, timeout=300)
        except Exception:  # noqa: BLE001
            return "fail"
        if len(raw) > 60e6:
            return "fail-too-large"
        open(setpath, "wb").write(raw)
    key = contig.split(".")[0]
    seq, take, hdr = [], False, ""
    for line in gzip.open(setpath, "rt"):
        if line.startswith(">"):
            if take:
                break
            take = key in line
            hdr = line
            continue
        if take:
            seq.append(line.strip())
    s = "".join(seq)
    if not s:
        return "fail"
    with gzip.open(path, "wt") as f:
        f.write(f">{contig}:{start}-{end} from ENA set {setname} contig_length={len(s)}\n{s[start - 1:end]}\n")
    return f"ena:{len(s)}"


todo = [r for r in stage0 if r["md5"] not in done][:MAX]
print("to resolve:", len(todo), flush=True)
out = open(SPECS, "a")
for r in todo:
    spec = {"md5": r["md5"], "rep_accession": r["rep_accession"], "length": int(r["length"]),
            "organism": r["organism"], "superkingdom": r["superkingdom"], "source": "", "contig": "",
            "strand": "", "segments": [], "n_contigs_in_model": 0, "note": ""}
    try:
        accs = [a for a in r["ncbi_accessions"].split(";") if a]
        res = None
        for a in accs:                                    # genomic coded_by first
            cb = CODED.get(a, "")
            if cb and not cb.startswith(("XM_", "NM_")):
                res = parse_location(cb)
                spec["source"] = f"coded_by:{a}"
                break
        if res is None:
            for a in accs:
                if a.startswith(("XP_", "NP_")) and rec.get(a, {}).get("gene_id"):
                    g = from_gene(rec[a]["gene_id"], a)
                    if g:
                        res = g[:4]
                        spec["source"] = f"datasets-gene:{rec[a]['gene_id']}"
                        if g[4]:
                            spec["note"] = f"transcript has {g[4]} bases not in genomic exons"
                    break
        if res is None:
            for a in accs:
                if a.startswith("WP_"):
                    res = from_ipg(a)
                    spec["source"] = f"ipg:{a}"
                    break
        if res is None and r["upi"]:
            for x in xrefs.get(r["upi"].split(";")[0], []):
                if x.get("active") and x["database"] in ("EMBLWGS", "EMBL", "EMBL_CON", "EMBL_TSA"):
                    res = from_ena(x["id"])
                    spec["source"] = f"ena:{x['id']}"
                    if res:
                        break
        if res is None:
            spec["note"] = "locus not resolved"
        else:
            spec["contig"], spec["strand"], spec["segments"], spec["n_contigs_in_model"] = res
    except Exception as e:  # noqa: BLE001
        spec["note"] = f"error: {e!r}"[:200]
    out.write(json.dumps(spec) + "\n")
    out.flush()
    done[r["md5"]] = spec
    print(len(done), r["rep_accession"], spec["source"], spec["contig"], len(spec["segments"]), spec["note"], flush=True)
out.close()

# ---- DNA per locus (one download per contig region; isoforms share it) ----------------
specs = [s for s in done.values() if s["segments"]]
clen_path = os.path.join(LD, "contig_lengths.json")
clen = json.load(open(clen_path)) if os.path.exists(clen_path) else {}
need = sorted({s["contig"] for s in specs if s["contig"] not in clen})
clen.update(contig_lengths(need))
json.dump(clen, open(clen_path, "w"))

regions = {}
for s in specs:
    lo = min(a for a, b in s["segments"])
    hi = max(b for a, b in s["segments"])
    key = s["contig"]
    regions.setdefault(key, []).append((lo, hi, s["md5"]))
fetched = 0
index_path = os.path.join(LD, "dna_index.json")
index = json.load(open(index_path)) if os.path.exists(index_path) else {}
for contig, rs in regions.items():
    rs.sort()
    # merge overlapping spans on the same contig (isoforms of one gene)
    merged = []
    for lo, hi, m in rs:
        if merged and lo <= merged[-1][1] + FLANK:
            merged[-1][1] = max(merged[-1][1], hi)
            merged[-1][2].append(m)
        else:
            merged.append([lo, hi, [m]])
    for lo, hi, ms in merged:
        L = clen.get(contig)
        a, b = max(1, lo - FLANK), (min(L, hi + FLANK) if L else hi + FLANK)
        lid = f"{contig}_{a}_{b}"
        path = os.path.join(LD, "dna", lid + ".fa.gz")
        if lid not in index or not os.path.exists(path):
            how = fetch_dna(contig, a, b, path)
            if how.startswith("ena:"):
                clen[contig] = int(how.split(":")[1])
            if contig not in clen and os.path.exists(path):
                hdr = gzip.open(path, "rt").readline()
                mm = re.search(r"length_(\d+)_cov", hdr)       # SPAdes contig name
                if mm:
                    clen[contig] = int(mm.group(1))
                    how += "+len-from-contig-name"
                else:                                           # last resort: full contig
                    try:
                        t = get(EU + "efetch.fcgi", {"db": "nuccore", "id": contig, "rettype": "fasta",
                                                     "retmode": "text"}, timeout=300).decode()
                        clen[contig] = len("".join(t.split("\n")[1:]).replace("\r", ""))
                        how += "+len-from-full-fetch"
                    except Exception:  # noqa: BLE001
                        pass
            index[lid] = {"contig": contig, "start": a, "end": b, "how": how, "md5s": ms}
            fetched += 1
            json.dump(index, open(index_path, "w"))
            json.dump(clen, open(clen_path, "w"))
            print("dna", lid, how, flush=True)
        else:
            index[lid]["md5s"] = sorted(set(index[lid]["md5s"]) | set(ms))
json.dump(index, open(index_path, "w"))
print("specs", len(done), "with segments", len(specs), "regions fetched now", fetched)
