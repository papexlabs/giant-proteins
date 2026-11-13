"""Positive and negative controls for the source-read check (cloud/dna_check.py), built from the finished
genome of Nocardia sp. NBC_01730 (CP109162.1, GCA_035920445.1; Illumina + ONT, complete chromosome).

The truth is known from outside our pipeline: the published complete genome is the reference truth and
the modifications below are introduced by us into a copy of it, so the real reads (SRR21428551 Illumina,
SRR21563524 ONT) do not support them. Regions were chosen from a self-alignment of the genome
(minimap2 -x asm5 -DP, repeats >= 2 kb at >= 98 % identity) so that each test region contains no other
long repeat.

  (a) deletion misjoin: 10 kb of unique sequence removed (genome 4,305,001-4,315,000); test region
      4,250,001-4,370,000 of the genome (110 kb after the deletion)
  (b) collapsed repeat: one of the two copies of a 4.55-kb 100 %-identical repeat removed
      (genome 7,542,232-7,546,784; the remaining copy is 89,952-94,504); test region 40,001-150,000
  (c) expanded repeat: a 5-kb unique segment (genome 5,858,001-5,863,000) tandem-duplicated, so the
      assembly carries an identical extra copy the genome lacks; test region 5,800,001-5,920,000
  (d) true negative: the unmodified genome around one of three rRNA-operon copies (5.5-5.7 kb,
      genome 263,060-268,735); test region 210,001-330,000
  (e) true negative: the unmodified genome around the 13.8-kb segmental repeat at the start of the
      NBC_01730 giant ORF (copies 2,943,113-2,956,904 and 6,933,167-6,946,964); test region
      6,880,001-7,000,000
Writes controls/dna/refs/{modified,original}.fa.gz, controls/dna/refs/modifications.json and the job
files controls/dna/jobs/*.json (Illumina only, and Illumina + ONT).
"""
import gzip
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
GEN = os.path.join(ROOT, "data", "loci", "nocardia_dl", "GCA_035920445.1.fna.gz")
READS = os.path.join(HERE, "reads")
ILL = [os.path.join(READS, "SRR21428551_1.fastq.gz"), os.path.join(READS, "SRR21428551_2.fastq.gz")]
ONT = [os.path.join(READS, "SRR21563524_1.fastq.gz")]
ILL_N, ONT_N = 15565314, 36579          # ENA read counts as used by the original NBC_01730 check
os.makedirs(os.path.join(HERE, "refs"), exist_ok=True)
os.makedirs(os.path.join(HERE, "jobs"), exist_ok=True)

name, buf = None, []
for line in gzip.open(GEN, "rt"):
    if line[0] == ">":
        name = line[1:].split()[0]
    else:
        buf.append(line.strip().upper())
g = "".join(buf)
assert name == "CP109162.1", name

# modifications, 1-based inclusive genome coordinates, applied from the highest position down so
# that each earlier (upstream) coordinate stays valid
DEL_A = (4305001, 4315000)
COLLAPSE_B = (7542232, 7546784)          # copy removed; kept copy 89,952-94,504
DUP_C = (5858001, 5863000)               # duplicated right after itself
assert g[89951:94504] == g[COLLAPSE_B[0] - 1:COLLAPSE_B[1]], "the two repeat copies are not identical"
m = g
m = m[:COLLAPSE_B[0] - 1] + m[COLLAPSE_B[1]:]
seg = m[DUP_C[0] - 1:DUP_C[1]]
m = m[:DUP_C[1]] + seg + m[DUP_C[1]:]
m = m[:DEL_A[0] - 1] + m[DEL_A[1]:]
la = DEL_A[1] - DEL_A[0] + 1


def mod_pos(p):
    """genome position (not inside a removed stretch) -> position in the modified sequence"""
    q = p
    if p > DEL_A[1]:
        q -= la
    if p > DUP_C[1]:
        q += DUP_C[1] - DUP_C[0] + 1
    return q


regions = {
    "a_deletion": {"ref": "modified", "genome_region": (4250001, 4370000),
                   "region": (4250001, 4370000 - la), "junction": DEL_A[0] - 1,
                   "truth": "misassembly: 10 kb of the genome missing (deletion misjoin at the junction)"},
    "b_collapse": {"ref": "modified", "genome_region": (40001, 150000), "region": (40001, 150000),
                   "repeat_copy": (89952, 94504),
                   "truth": "misassembly: a two-copy 4.55-kb repeat collapsed to one copy (reads of both copies pile here)"},
    "c_expansion": {"ref": "modified", "genome_region": (5800001, 5920000),
                    "region": (mod_pos(5800001), mod_pos(5920000)),
                    "copies": ((DUP_C[0] - la, DUP_C[1] - la),
                               (DUP_C[1] - la + 1, DUP_C[1] - la + (DUP_C[1] - DUP_C[0] + 1))),
                    "truth": "misassembly: a 5-kb unique segment present twice in tandem (extra copy not in the genome)"},
    "d_rrna_correct": {"ref": "original", "genome_region": (210001, 330000), "region": (210001, 330000),
                       "repeat_copy": (263060, 268735),
                       "truth": "correct assembly; contains one copy of a three-copy rRNA operon repeat"},
    "e_segdup_correct": {"ref": "original", "genome_region": (6880001, 7000000), "region": (6880001, 7000000),
                         "repeat_copy": (6933167, 6946964),
                         "truth": "correct assembly; contains a 13.8-kb segmental repeat (second copy at 2.94 Mb, reversed)"},
}
# sanity: the modified sequence around the deletion junction and the duplication
assert m[mod_pos(DEL_A[0] - 100) - 1:mod_pos(DEL_A[0] - 100) + 99] == g[DEL_A[0] - 101:DEL_A[0] - 1]
c1, c2 = regions["c_expansion"]["copies"]
assert m[c1[0] - 1:c1[1]] == m[c2[0] - 1:c2[1]] == seg, "duplication copies"

for tag, seq in (("modified", m), ("original", g)):
    p = os.path.join(HERE, "refs", f"{tag}.fa.gz")
    if not os.path.exists(p):
        with gzip.open(p, "wt", compresslevel=3) as f:
            f.write(">CP109162.1\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
json.dump({"genome": "GCA_035920445.1 CP109162.1", "genome_len": len(g), "modified_len": len(m),
           "deletion_a": DEL_A, "collapse_b_removed_copy": COLLAPSE_B, "duplication_c": DUP_C, "regions": regions},
          open(os.path.join(HERE, "refs", "modifications.json"), "w"), indent=1)

for tag, r in regions.items():
    for mode in ("ill", "illont"):
        # paths relative to the repository root (run_controls.sh runs from there): no absolute home paths in jobs
        rel = lambda x: os.path.relpath(x, ROOT)  # noqa: E731
        job = {"job_id": f"ctl_{tag}_{mode}", "assembly_path": rel(os.path.join(HERE, "refs", r["ref"] + ".fa.gz")),
               "contig": "CP109162.1", "gene_start": r["region"][0], "gene_end": r["region"][1], "strand": "+",
               "illumina": [rel(x) for x in ILL], "illumina_expected_reads": ILL_N,
               "ont": [rel(x) for x in ONT] if mode == "illont" else [],
               "ont_expected_reads": ONT_N, "truth": r["truth"]}
        json.dump(job, open(os.path.join(HERE, "jobs", job["job_id"] + ".json"), "w"), indent=1)
print(json.dumps({k: v["region"] for k, v in regions.items()}))
print("genome", len(g), "modified", len(m))
