"""Every number quoted in the manuscript, computed from the repository's own files.

Writes paper/numbers.json (key -> value), paper/NUMBERS_generated.tsv (key, value, source) and
paper/numbers.tex (one \\nmdef{key}{value} per key; the LaTeX sources print a number only through
\\nm{key}, so every printed number is one of these values). paper/NUMBERS.md maps the text to the keys.
Run from anywhere: python3 paper/make_numbers.py   (standard library only)
"""
import csv
import glob
import gzip
import io
import json
import os
import re
import statistics
import tarfile
from collections import Counter

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
D = os.path.join(ROOT, "data")
N = {}          # key -> (value, source)


def hu(x, nd=0):
    """Round half up (Python's round() rounds 72.5 to 72)."""
    import decimal
    q = decimal.Decimal(1).scaleb(-nd)
    v = decimal.Decimal(repr(x)).quantize(q, rounding=decimal.ROUND_HALF_UP)
    return int(v) if nd == 0 else float(v)


def put(key, value, source):
    N[key] = (value, source)


def tsv(name):
    return list(csv.DictReader(open(os.path.join(D, name)), delimiter="\t"))


TITIN, PKZ = 35991, 45212
put("titin_len", TITIN, "human titin NP_001254479.2 length, data/locus_checks.tsv (control row, column length)")
put("pkz_len", PKZ, "PKZILLA-1 KAL1508129.1 length, data/locus_checks.tsv column length")

# ---- inventory ---------------------------------------------------------------------------
up = tsv("uniprotkb_ge35992.tsv")
put("uniprotkb_gt_titin", len(up), "data/uniprotkb_ge35992.tsv rows (UniProtKB 2026_04, length > 35,991)")
put("uniprotkb_gt_pkz", sum(1 for r in up if int(r["Length"]) > PKZ), "data/uniprotkb_ge35992.tsv, Length > 45,212")
pe = Counter(r["Protein existence"] for r in up)
put("uniprotkb_pe3", pe["Inferred from homology"], "data/uniprotkb_ge35992.tsv column 'Protein existence' = Inferred from homology")
put("uniprotkb_pe4", pe["Predicted"], "data/uniprotkb_ge35992.tsv column 'Protein existence' = Predicted")
put("uniprotkb_reviewed", sum(1 for r in up if r["Reviewed"] == "reviewed"), "data/uniprotkb_ge35992.tsv column Reviewed")
upa = tsv("uniparc_ge35992.tsv")
put("uniparc_gt_titin", len(upa), "data/uniparc_ge35992.tsv rows")
put("uniparc_gt_pkz", sum(1 for r in upa if int(r["length"]) > PKZ), "data/uniparc_ge35992.tsv, length > 45,212")
nc = {r["query_name"]: int(r["count"]) for r in tsv("ncbi_counts.tsv")}
put("ncbi_gt_titin", nc["all_gt_titin_isoform"], "data/ncbi_counts.tsv all_gt_titin_isoform (Entrez, 11 Nov 2026)")
put("ncbi_gt_pkz", nc["all_gt_pkzilla1"], "data/ncbi_counts.tsv all_gt_pkzilla1")
put("refseq_gt_titin", nc["refseq_gt_titin_isoform"], "data/ncbi_counts.tsv refseq_gt_titin_isoform")
g = tsv("giant_unique_sequences.tsv")
put("union_seqs", len(g), "data/giant_unique_sequences.tsv rows (unique MD5)")
put("union_gt_pkz", sum(1 for r in g if int(r["length"]) > PKZ), "data/giant_unique_sequences.tsv, length > 45,212")
put("union_live", sum(1 for r in g if r["currently_live"] == "yes"), "data/giant_unique_sequences.tsv currently_live = yes")
put("excl_nolive", sum(1 for r in g if r["currently_live"] == "no"), "data/giant_unique_sequences.tsv currently_live = no")
put("excl_patent", sum(1 for r in g if r["patent_or_synthetic"] == "yes"), "data/giant_unique_sequences.tsv patent_or_synthetic = yes")
put("excl_partial", sum(1 for r in g if r["partial_in_title"] == "yes"), "data/giant_unique_sequences.tsv partial_in_title = yes")
put("longest_record", max(int(r["length"]) for r in g), "data/giant_unique_sequences.tsv max length (a patent sequence, excluded)")
s0 = tsv("stage0_table.tsv")
kept = [r for r in s0 if r["triage"] != "T-excl"]
put("n_retained", len(kept), "data/stage0_table.tsv triage != T-excl")
put("n_excluded", len(s0) - len(kept), "data/stage0_table.tsv triage = T-excl")
put("n_retained_gt_pkz", sum(1 for r in kept if int(r["length"]) > PKZ), "data/stage0_table.tsv, retained and length > 45,212")
loci0 = tsv("stage0_loci.tsv")
put("loci_retained_stage0", sum(1 for r in loci0 if r["triage"] != "T-excl"), "data/stage0_loci.tsv triage != T-excl")
xp = [r for r in tsv("ncbi_records_parsed.tsv") if r["accession"].startswith("XP_")]
put("xp_records", len(xp), "data/ncbi_records_parsed.tsv accessions starting XP_")
put("xp_genes", len({r["gene_id"] for r in xp if r["gene_id"]}), "data/ncbi_records_parsed.tsv distinct gene_id among XP_")

# ---- grades ------------------------------------------------------------------------------
lc = tsv("locus_checks.tsv")
main = [r for r in lc if not r["md5"].startswith("control_")]
big = [r for r in main if int(r["length"]) > PKZ]
sk = Counter(r["superkingdom"] for r in main)
put("n_euk", sk["Eukaryota"], "data/locus_checks.tsv superkingdom (510 rows, controls excluded)")
put("n_bact", sk["Bacteria"], "data/locus_checks.tsv superkingdom")
put("n_arch", sk["Archaea"], "data/locus_checks.tsv superkingdom")
skb = Counter(r["superkingdom"] for r in big)
put("big_euk", skb["Eukaryota"], "data/locus_checks.tsv superkingdom, length > 45,212")
put("big_bact", skb["Bacteria"], "data/locus_checks.tsv superkingdom, length > 45,212")
for col, tag in (("grade_stage1", "s1"), ("grade_stage2_rna", "s2"), ("grade", "fin")):
    for scope, rows in (("all", main), ("big", big)):
        c = Counter(r[col] for r in rows)
        for gr in "ABCD":
            put(f"{tag}_{scope}_{gr}", c[gr], f"data/locus_checks.tsv column {col}, {'510 sequences' if scope == 'all' else 'length > 45,212'}")
tr = Counter((r["grade_stage1"], r["grade_stage2_rna"]) for r in main)
put("s1s2_C_to_B", tr[("C", "B")], "data/locus_checks.tsv grade_stage1 C -> grade_stage2_rna B")
put("s1s2_B_to_A", tr[("B", "A")], "data/locus_checks.tsv grade_stage1 B -> grade_stage2_rna A")
apo = sum(1 for r in main if r["grade_stage1"] == "C" and r["grade_stage2_rna"] == "B" and r["organism"] == "Apostichopus japonicus")
put("s1s2_C_to_B_apostichopus", apo, "data/locus_checks.tsv, C -> B rows with organism Apostichopus japonicus")
dr, drb = Counter(), Counter()
for r in main:
    if r["grade_stage1"] == "D":
        for k in ("D1", "D2", "D3", "D4", "D5"):
            if (k + " repeat translation" if k == "D1" else k + " ") in r["grade_reason_stage1"]:
                dr[k] += 1
                if int(r["length"]) > PKZ:
                    drb[k] += 1
for k in ("D1", "D2", "D3", "D5"):
    put(f"rule_{k}_all", dr[k], f"data/locus_checks.tsv grade_reason_stage1 contains '{k} ' (D rows)")
    put(f"rule_{k}_big", drb[k], f"same, length > 45,212")
d1_euk_big = sum(1 for r in big if "D1 repeat translation" in r["grade_reason_stage1"] and r["superkingdom"] == "Eukaryota")
put("rule_D1_big_euk", d1_euk_big, "data/locus_checks.tsv D1 rows > 45,212 aa, superkingdom Eukaryota")

# the D1 pattern (repeat translation) with and without a DNA-level contradiction (D1 rule v2)
pat = [r for r in main if "D1 repeat translation" in r["grade_reason_stage1"] or "repeat flag" in r["grade_reason_stage1"]]
put("d1_pattern_all", len(pat), "data/locus_checks.tsv grade_reason_stage1 contains 'D1 ' or 'repeat flag'")
put("d1_pattern_big", sum(1 for r in pat if int(r["length"]) > PKZ), "same, length > 45,212")
put("d1_pattern_big_euk", sum(1 for r in pat if int(r["length"]) > PKZ and r["superkingdom"] == "Eukaryota"), "same, Eukaryota")
put("repeat_flag_all", sum(1 for r in main if "repeat flag" in r["grade_reason_stage1"]), "data/locus_checks.tsv 'repeat flag' (D1 pattern without D2/D3)")
put("repeat_flag_big", sum(1 for r in big if "repeat flag" in r["grade_reason_stage1"]), "same, length > 45,212")
fb = [r for r in big if "repeat flag" in r["grade_reason_stage1"]]
weak = [float(r["fickett_median"]) for r in fb if r["fickett_median"] and float(r["fickett_median"]) < 0.95]
put("flag_big_weak_coding", len(weak), "data/locus_checks.tsv repeat-flag rows > 45,212 aa with Fickett median < 0.95")
put("flag_big_weak_coding_range", (min(weak), max(weak)), "same, Fickett median range")
put("flag_big_coding", sum(1 for r in fb if r["fickett_median"] and float(r["fickett_median"]) >= 0.95), "repeat-flag rows > 45,212 aa with Fickett median >= 0.95")
put("flag_big_nodna", sum(1 for r in fb if not r["fickett_median"]), "repeat-flag rows > 45,212 aa without locus DNA")
put("d_big_d1_with_dna", sum(1 for r in big if r["grade"] == "D" and "D1 repeat translation" in r["grade_reason_stage1"]), "data/locus_checks.tsv D rows > 45,212 aa with D1 (now always with D2 or D3)")
put("d_big_d5", sum(1 for r in big if r["grade"] == "D" and "D5 " in r["grade_reason_stage1"]), "data/locus_checks.tsv D rows > 45,212 aa with D5")
put("rna_not_tested_d1", sum(1 for r in main if r["rnaseq_verdict"].startswith("not-tested (graded D1")), "data/locus_checks.tsv rnaseq_verdict not-tested (graded D1 ...)")

# ---- RNA-seq stage of the repeat-flag models (cloud/bundle_rep_rna, 13 Nov) -------------------------
REP_LOCI = {}
for pth in glob.glob(os.path.join(ROOT, "cloud", "bundle_rep_rna", "loci", "*.json")):
    for L in json.load(open(pth)):
        REP_LOCI[L["md5"]] = os.path.basename(pth)[:-5]
rep_jobs = list(csv.DictReader(open(os.path.join(ROOT, "cloud", "bundle_rep_rna", "jobs.tsv")), delimiter="\t"))
rep_job_sp = {j["species_key"] for j in rep_jobs}
rep_rows = [r for r in main if r["md5"] in REP_LOCI or ("repeat flag" in r["grade_reason_stage1"] and int(r["length"]) and
            r["rnaseq_verdict"].startswith("not-tested (locus not resolved"))]
put("rep_models", len([r for r in main if "repeat flag" in r["grade_reason_stage1"] and
                       (r["md5"] in REP_LOCI or r["rnaseq_verdict"].startswith("not-tested (locus not resolved"))]),
    "data/locus_checks.tsv repeat-flag models in cloud/bundle_rep_rna/loci or without a locus (the 43 never sent to the RNA stage)")
tested = [r for r in main if REP_LOCI.get(r["md5"]) in rep_job_sp]
put("rep_testable", len(tested), "data/locus_checks.tsv models whose bundle_rep_rna species key has runs in jobs.tsv")
put("rep_testable_big", sum(1 for r in tested if int(r["length"]) > PKZ), "same, length > 45,212")
put("rep_jobs", len(rep_jobs), "cloud/bundle_rep_rna/jobs.tsv rows")
put("rep_runs", sum(int(j["n_runs"]) for j in rep_jobs), "cloud/bundle_rep_rna/jobs.tsv n_runs sum")
put("rep_gb", round(sum(float(j["runs_gb"]) for j in rep_jobs)), "cloud/bundle_rep_rna/jobs.tsv runs_gb sum")
rsup = {r["md5"]: r for r in tsv("rnaseq_support.tsv")}
for acc, key in (("CAH8664242.1", "schisto1"), ("CAH8664230.1", "schisto2"), ("CAH1773861.1", "owenia"), ("CAR9444729.1", "arcto")):
    r = rsup[next(m for m in rsup if rsup[m]["rep_accession"] == acc)]
    put(f"rep_{key}_supported", int(r["introns_counted_supported"]), f"data/rnaseq_support.tsv {acc} introns_counted_supported")
    put(f"rep_{key}_introns", int(r["n_introns"]), f"data/rnaseq_support.tsv {acc} n_introns")
    put(f"rep_{key}_testable", int(r["testable_introns"]), f"data/rnaseq_support.tsv {acc} testable_introns")
    put(f"rep_{key}_testable_unsup", int(r["testable_introns_unsupported"]), f"data/rnaseq_support.tsv {acc} testable_introns_unsupported")
    put(f"rep_{key}_unique_reads", int(r["reads_unique"]), f"data/rnaseq_support.tsv {acc} reads_unique")
    put(f"rep_{key}_verdict", r["verdict"], f"data/rnaseq_support.tsv {acc} verdict")
    put(f"rep_{key}_grade", next(x["grade"] for x in main if x["rep_accession"] == acc), f"data/locus_checks.tsv {acc} grade")
noexp = [r for r in tested if rsup.get(r["md5"], {}).get("verdict") == "no-expression"]
put("rep_noexpr", len(noexp), "data/rnaseq_support.tsv verdict no-expression among the testable repeat-flag models")
put("rep_noexpr_max_unique", max(int(rsup[r["md5"]]["reads_unique"]) for r in noexp), "same, maximum reads_unique")
put("rep_contradicted", sum(1 for r in tested if rsup.get(r["md5"], {}).get("verdict") == "contradicted"), "same, verdict contradicted")
put("rep_to_B", sum(1 for r in tested if r["grade"] == "B" and r["grade_stage1"] == "C"), "data/locus_checks.tsv testable repeat-flag models C -> B")
put("rep_to_D", sum(1 for r in tested if r["grade"] == "D"), "same, grade D")
# run results of that stage, read from the stored archive (the synced cloud_out/ is git-ignored)
with tarfile.open(os.path.join(D, "rnaseq_run_results.tar.gz")) as _t:
    _rep = {m.name: m for m in _t.getmembers() if m.isfile() and "_rep/" in m.name and m.name.endswith(".json")}
    put("rep_run_results", sum(1 for n in _rep if not n.endswith(".skip.json")),
        "data/rnaseq_run_results.tar.gz results/*_rep/*.json (run results, skips excluded)")
    put("rep_run_skips", sum(1 for n in _rep if n.endswith(".skip.json")), "data/rnaseq_run_results.tar.gz results/*_rep/*.skip.json")
    seen = json.load(_t.extractfile(_rep["results/Schistosoma_haematobium_944470465.2_rep/ERR10016938.json"]))
put("rep_schisto_err10016938", next(L["introns_with_ge3_exact_reads"] for L in seen["loci"] if L["rep_accession"] == "CAH8664242.1"),
    "data/rnaseq_run_results.tar.gz results/Schistosoma_haematobium_944470465.2_rep/ERR10016938.json introns_with_ge3_exact_reads "
    "(CAH8664242.1), a run not seen during the Mac test")
mp = json.load(open(os.path.join(ROOT, "controls", "rna_rep", "junction_mappability.json")))
mpb = {x["rep_accession"]: x for x in mp if x.get("n_introns")}
for acc, key in (("CAR0409254.1", "gas1"), ("CAR0409140.1", "gas2"), ("CAR0409152.1", "gas3"), ("KAL7845143.1", "acn"),
                 ("KAJ8332361.1", "syn1"), ("KAJ8332362.1", "syn2"), ("CAH8664242.1", "schisto1"), ("CAH1773861.1", "owenia"),
                 ("CAR9444729.1", "arcto")):
    put(f"mapp_{key}_uj", mpb[acc]["introns_unique_junction"], f"controls/rna_rep/junction_mappability.json {acc} introns_unique_junction")
    put(f"mapp_{key}_n", mpb[acc]["n_introns"], f"controls/rna_rep/junction_mappability.json {acc} n_introns")
    put(f"mapp_{key}_pct", round(100 * mpb[acc]["share_unique_junction"]), f"same, share_unique_junction x 100")

# ---- controls ------------------------------------------------------------------------------
CT = os.path.join(ROOT, "controls")
c1 = json.load(open(os.path.join(CT, "d1", "summary.json")))
put("ctl_release", c1["uniprot_release"], "controls/d1/summary.json uniprot_release")
put("ctl_n", c1["n_entries"], "controls/d1/summary.json n_entries")
for g, key in (("G1 silk (Swiss-Prot)", "silk"), ("G1b complete spidroins", "spidroin"), ("G2 mucins", "mucin"),
               ("G3 bacterial adhesins", "adhesin"), ("G4 long proteins", "long")):
    put(f"ctl_{key}_n", c1["by_group"][g]["n"], f"controls/d1/summary.json by_group[{g}].n")
    put(f"ctl_{key}_d1", c1["by_group"][g]["D1"], f"controls/d1/summary.json by_group[{g}].D1")
put("ctl_d1_flagged", c1["all"]["D1"], "controls/d1/summary.json all.D1")
put("ctl_d23_flagged", sum(1 for x in c1["all"]["flagged"] if x[8] or x[9]), "controls/d1/summary.json flagged with D2 or D3")
ctl_rows = list(csv.DictReader(open(os.path.join(CT, "d1", "d1_controls.tsv")), delimiter="\t"))
fl = [r for r in ctl_rows if r["D1"] == "True"]
put("ctl_flagged_with_cds", sum(1 for r in fl if r.get("cds_protein_id")), "controls/d1/d1_controls.tsv D1 rows with a matched CDS")
put("ctl_flagged_fickett_min", min(float(r["fickett_median"]) for r in fl if r.get("fickett_median")), "same, minimum Fickett median")
put("ctl_flagged_open_max", max(int(r["n_other_frames_open"]) for r in fl if r.get("n_other_frames_open")), "same, maximum open alternative frames")
hp = os.path.join(CT, "d1", "holdout_summary.json")
if os.path.exists(hp):
    h = json.load(open(hp))
    put("hold_n", h["n_holdout"], "controls/d1/holdout_summary.json n_holdout")
    put("hold_trigger", h["repeat_trigger"], "controls/d1/holdout_summary.json repeat_trigger")
    put("hold_d1_v1", h["D1_rule_v1"], "controls/d1/holdout_summary.json D1_rule_v1")
    put("hold_d1_cds", h["D1_with_cds"], "controls/d1/holdout_summary.json D1_with_cds")
    put("hold_d1_v2", h["D1_revised_D1_and_D2_or_D3"], "controls/d1/holdout_summary.json D1_revised")
pw = json.load(open(os.path.join(CT, "rna_power", "summary.json")))
put("power_models", pw["models"], "controls/rna_power/summary.json models")
put("power_reps", pw["replicates"], "controls/rna_power/summary.json replicates")
for k in ("shift 0.1", "shift 0.2", "shift 0.3", "shift 0.4", "shift 0.5", "shift 0.6", "skip 0.2", "skip 0.3", "skip 0.5", "skip 0.6"):
    for o in ("v1", "v2"):
        for m in ("contradicted", "under-supported where covered", "partial", "full"):
            put(f"power_{k.replace(' ', '_')}_{o}_{m.split()[0]}", pw["by_kind_fraction"][k][o][m],
                f"controls/rna_power/summary.json by_kind_fraction['{k}'].{o}['{m}']")
bk = pw["by_kind_fraction"]


def down(k, o):
    return bk[k][o]["contradicted"] + bk[k][o]["under-supported where covered"]


SPW = "controls/rna_power/summary.json by_kind_fraction"
put("pw_contr_ge60_min_pct", hu(100 * min(bk[f"{t} {f}"]["v2"]["contradicted"] for t in ("shift", "skip") for f in ("0.6", "0.8"))),
    SPW + ": minimum 'contradicted' share at 60 and 80 % wrong introns (shift, skip)")
put("pw_shift50_pct", hu(100 * bk["shift 0.5"]["v2"]["contradicted"]), SPW + "['shift 0.5'].v2.contradicted")
put("pw_skip50_pct", hu(100 * bk["skip 0.5"]["v2"]["contradicted"]), SPW + "['skip 0.5'].v2.contradicted")
put("pw_v1_partial_2040_min_pct", hu(100 * min(bk[f"{t} {f}"]["v1"]["partial"] for t in ("shift", "skip") for f in ("0.2", "0.3", "0.4")), 1),
    SPW + ": minimum v1 'partial' share at 20-40 % wrong")
put("pw_v2_down_3040_min_pct", hu(100 * min(down(f"{t} {f}", "v2") for t in ("shift", "skip") for f in ("0.3", "0.4"))),
    SPW + ": minimum v2 downgraded share at 30-40 % wrong")
put("pw_v2_shift20_pct", hu(100 * down("shift 0.2", "v2")), SPW + "['shift 0.2'].v2 contradicted + under-supported")
put("pw_v2_skip20_pct", hu(100 * down("skip 0.2", "v2")), SPW + "['skip 0.2'].v2 contradicted + under-supported")
put("pw_v2_10_down_pct", hu(100 * max(down(f"{t} 0.1", "v2") for t in ("shift", "skip"))), SPW + ": maximum v2 downgraded share at 10 % wrong")
dn = json.load(open(os.path.join(CT, "dna", "summary.json")))
for k, v in dn["v3"].items():
    put(f"dnactl_v3_{k.replace('/', '_')}", v["final_verdict"], f"controls/dna/summary.json v3.{k}.final_verdict")
for k, v in dn["v3.1"].items():
    put(f"dnactl_v31_{k.replace('/', '_')}", v["final_verdict"], f"controls/dna/summary.json v3.1.{k}.final_verdict")
    if v.get("copy_number_blocks"):
        put(f"dnactl_v31_{k.replace('/', '_')}_ratio", v["copy_number_blocks"][0]["depth_ratio"], f"controls/dna/summary.json v3.1.{k}.copy_number_blocks")
put("dnactl_real_unchanged", all(v["verdict_v3"] == v["verdict_v3.1"] for v in dn["rejudged_real_v3.1"].values()),
    "controls/dna/summary.json rejudged_real_v3.1 (v3 verdict == v3.1 verdict for all five)")
rr = [x[2] for v in dn["rejudged_real_v3.1"].values() for x in v.get("repeat_block_all_mapq_depth_ratio", [])]
put("dnactl_repeat_depth_ratio_range", (min(rr), max(rr)), "controls/dna/summary.json repeat_block_all_mapq_depth_ratio (all five checks)")
omni = [x[2] for k, v in dn["rejudged_real_v3.1"].items() if k.startswith("omni") for x in v.get("repeat_block_all_mapq_depth_ratio", [])]
put("dnactl_omni_repeat_depth_ratio_range", (min(omni), max(omni)), "same, Omnitrophota only")
sens = list(csv.DictReader(open(os.path.join(CT, "sensitivity", "sensitivity_stage12.tsv")), delimiter="\t"))
put("sens_variants", len(sens) - 1, "controls/sensitivity/sensitivity_stage12.tsv rows minus the baseline")
sbase = sens[0]
put("sens_big_ab_stable", sum(1 for r in sens[1:] if r["big_A"] == sbase["big_A"] and r["big_B"] == sbase["big_B"]),
    "controls/sensitivity/sensitivity_stage12.tsv variant rows with the baseline A and B counts above 45,212 aa")
put("sens_big_ab_changed", [(r["variant"], r["big_A"], r["big_B"]) for r in sens[1:] if (r["big_A"], r["big_B"]) != (sbase["big_A"], sbase["big_B"])],
    "controls/sensitivity/sensitivity_stage12.tsv variants whose A or B count above 45,212 aa differs from the baseline")
s6 = next(r for r in sens if r["variant"].startswith("S6 peptides"))
put("sens_s6_lenient_all_B", int(s6["all_B"]), "controls/sensitivity/sensitivity_stage12.tsv S6 >= 3 fifths variant, all_B")
put("sens_base_all_B", int(sbase["all_B"]), "controls/sensitivity/sensitivity_stage12.tsv baseline all_B")
d1v1 = next(r for r in sens if r["variant"].startswith("D1 rule v1:"))
put("sens_d1v1_big_D", int(d1v1["big_D"]), "controls/sensitivity/sensitivity_stage12.tsv D1 rule v1 big_D")
put("sens_d1v1_all_D", int(d1v1["all_D"]), "same, all_D")
put("sens_allA_range", (min(int(r["all_A"]) for r in sens), max(int(r["all_A"]) for r in sens)), "controls/sensitivity/sensitivity_stage12.tsv all_A range")
dsens = list(csv.DictReader(open(os.path.join(CT, "dna", "sensitivity.tsv")), delimiter="\t"))
put("dna_sens_settings", len({r["setting"] for r in dsens}) - 1, "controls/dna/sensitivity.tsv settings minus default")
put("dna_sens_grade_changes", [(r["setting"], r["accession"], r["grade_now"], r["grade"]) for r in dsens if r["grade"] != r["grade_now"]],
    "controls/dna/sensitivity.tsv rows whose grade differs from the current grade")
# ---- proteomics (PRIDE PXD025008 re-search) ---------------------------------------------------
PR = os.path.join(ROOT, "proteomics", "results")
p6 = json.load(open(os.path.join(PR, "full6_summary.json")))
pv = json.load(open(os.path.join(PR, "verify_extra.json")))
pf = json.load(open(os.path.join(PR, "summary.json")))
S6 = "proteomics/results/full6_summary.json"
SV = "proteomics/results/verify_extra.json"
put("prot_spectra", pv["spectra_in_peak_lists"], SV + " spectra_in_peak_lists (BEGIN IONS in the six MGF files)")
put("prot_db_entries", pv["database_entries"], SV + " database_entries")
put("prot_fractions", pv["fractions"], SV + " fractions (MudPIT fractions in the six peak lists)")
put("prot_scan_share_pct", round(100 * pv["share_of_scan_numbers_present"]), SV + " share_of_scan_numbers_present x 100")
put("prot_id_rate_pct", round(100 * p6["target_psms_q01"] / pv["spectra_in_peak_lists"], 1), S6 + " target_psms_q01 / " + SV + " spectra_in_peak_lists")
put("prot_distal", [(x[1], x[2]) for x in pv["distal_peptides_psms"]], SV + " distal_peptides_psms (start, PSMs) for starts > 6,000")
put("prot_psms", p6["target_psms_q01"], S6 + " target_psms_q01 (spectrum q <= 0.01)")
put("prot_decoy_psms", p6["decoy_psms_q01"], S6 + " decoy_psms_q01")
put("prot_peptides", p6["target_peptides_q01"], S6 + " target_peptides_q01 (peptide q <= 0.01)")
put("prot_proteins", p6["target_proteins_q01"], S6 + " target_proteins_q01 (protein q <= 0.01)")
gi = p6["giant"]
put("prot_giant_len", gi["length"], S6 + " giant.length")
put("prot_giant_peptides", gi["unique_peptides_q01"], S6 + " giant.unique_peptides_q01 (PSM, peptide and protein q <= 0.01; unique in the database)")
put("prot_giant_shared", gi["shared_peptides_q01"], S6 + " giant.shared_peptides_q01")
put("prot_giant_psms", gi["unique_psms_q01"], S6 + " giant.unique_psms_q01")
put("prot_giant_residues", gi["residues_covered_by_unique_peptides"], S6 + " giant.residues_covered_by_unique_peptides")
put("prot_giant_cov_pct", round(100 * gi["coverage_frac"], 1), S6 + " giant.coverage_frac x 100")
put("prot_giant_regions", gi["covered_regions"], S6 + " giant.covered_regions")
put("prot_giant_first", gi["first_unique_position"], S6 + " giant.first_unique_position")
put("prot_giant_last_end", pv["last_residue_covered"], SV + " last_residue_covered")
put("prot_giant_last_pct", round(100 * pv["fraction_of_length_up_to_last_residue"], 1), SV + " fraction_of_length_up_to_last_residue x 100")
put("prot_giant_fifths", gi["unique_peptides_per_fifth_of_length"], S6 + " giant.unique_peptides_per_fifth_of_length")
put("prot_giant_first6000", pv["unique_starts_first_6000"], SV + " unique_starts_first_6000")
put("prot_giant_after6000", pv["unique_starts_after_6000"], SV + " unique_starts_after_6000 (start positions)")
put("prot_giant_protein_q", pv["best_protein_q"]["giant"], SV + " best_protein_q.giant")
put("prot_rev_q", pv["best_protein_q"]["CTRLREV_giant"], SV + " best_protein_q.CTRLREV_giant")
put("prot_shuf_q", pv["best_protein_q"]["CTRLSHUF_giant"], SV + " best_protein_q.CTRLSHUF_giant")
put("prot_rev_psms", p6["controls_giant_copies"]["reversed giant"]["psms_q01"], S6 + " controls_giant_copies['reversed giant'].psms_q01")
put("prot_shuf_psms", p6["controls_giant_copies"]["shuffled giant"]["psms_q01"], S6 + " controls_giant_copies['shuffled giant'].psms_q01")
put("prot_tryptic_range", pv["theoretical_min_max_full_bins"], SV + " theoretical_min_max_full_bins (7-40 aa fully tryptic per full 1,000-aa bin)")
put("prot_repeat_span", pv["repeat_span"], SV + " repeat_span (20-mers occurring more than once, blocks >= 100 aa)")
put("prot_repeat_blocks", len(pv["repeat_blocks_ge100aa"]), SV + " repeat_blocks_ge100aa (count)")
put("prot_peptides_in_or_after_repeats", pv["unique_peptides_inside_or_after_first_block"], SV + " unique_peptides_inside_or_after_first_block")
for key, lab in (("eftu", "elongation factor Tu"), ("l7", "ribosomal protein L7/L12")):
    v = p6["normal_proteins"][lab]
    put(f"prot_{key}_peptides", v["unique_peptides_q01"], S6 + f" normal_proteins['{lab}'].unique_peptides_q01")
    put(f"prot_{key}_cov_pct", round(100 * v["coverage_by_unique_frac"]), S6 + f" normal_proteins['{lab}'].coverage_by_unique_frac x 100")
    put(f"prot_{key}_len", v["length"], S6 + f" normal_proteins['{lab}'].length")
ab = pf["abundant_proteins_same_organism"]
cov_ab = [x["unique_coverage_pct"] for x in ab if "WP_128698900.1" not in x["protein"]]
put("prot_top_n", len(ab), "proteomics/results/summary.json abundant_proteins_same_organism (the 12 Velamenicoccus proteins with most PSMs)")
put("prot_top_rank_giant", 1 + [x["protein"] for x in ab].index(next(x["protein"] for x in ab if "WP_128698900.1" in x["protein"])),
    "proteomics/results/summary.json abundant_proteins_same_organism, rank of the giant by PSMs")
put("prot_top_cov_range", (min(cov_ab), max(cov_ab)), "proteomics/results/summary.json abundant_proteins_same_organism unique_coverage_pct, the other 11")
put("prot_wall_min", round(pf["search_settings"]["wall_s"] / 60), "proteomics/results/summary.json search_settings.wall_s / 60")
put("prot_rss_gb", pf["search_settings"]["peak_rss_gb"], "proteomics/results/summary.json search_settings.peak_rss_gb")
# independent search of the 60 RAW files (13 Nov)
r6 = json.load(open(os.path.join(PR, "raw60_summary.json")))
vr = json.load(open(os.path.join(PR, "verify_raw.json")))
SR, SVR = "proteomics/results/raw60_summary.json", "proteomics/results/verify_raw.json"
put("raw_files", vr["raw_files"], SVR + " raw_files (converted, searched)")
put("raw_conversion_errors", vr["conversion_errors"], SVR + " conversion_errors")
put("raw_ms2", vr["ms2_spectra_raw"], SVR + " ms2_spectra_raw (ThermoRawFileParser metadata)")
put("raw_mgf_share_pct", round(100 * vr["mgf_share_of_raw_ms2"]), SVR + " mgf_share_of_raw_ms2 x 100")
put("raw_instrument", vr["instrument"][0], SVR + " instrument")
put("raw_psms", r6["target_psms_q01"], SR + " target_psms_q01")
put("raw_peptides", r6["target_peptides_q01"], SR + " target_peptides_q01")
put("raw_proteins", r6["target_proteins_q01"], SR + " target_proteins_q01")
rg = r6["giant"]
put("raw_giant_peptides", rg["unique_peptides_q01"], SR + " giant.unique_peptides_q01")
put("raw_giant_psms", rg["unique_psms_q01"], SR + " giant.unique_psms_q01")
put("raw_giant_residues", rg["residues_covered_by_unique_peptides"], SR + " giant.residues_covered_by_unique_peptides")
put("raw_giant_cov_pct", round(100 * rg["coverage_frac"], 1), SR + " giant.coverage_frac x 100")
put("raw_giant_fifths", rg["unique_peptides_per_fifth_of_length"], SR + " giant.unique_peptides_per_fifth_of_length")
put("raw_giant_first6000", vr["raw_starts_first_6000"], SVR + " raw_starts_first_6000")
put("raw_giant_after6000", [(x[0], x[1]) for x in vr["raw_after_6000"]], SVR + " raw_after_6000 (start, PSMs)")
put("raw_giant_last_end", vr["raw_last_residue"], SVR + " raw_last_residue")
put("raw_giant_last_pct", vr["raw_last_residue_pct"], SVR + " raw_last_residue_pct")
put("raw_giant_protein_q", vr["best_protein_q"]["giant"], SVR + " best_protein_q.giant")
put("raw_rev_q", vr["best_protein_q"]["CTRLREV_giant"], SVR + " best_protein_q.CTRLREV_giant")
put("raw_shuf_q", vr["best_protein_q"]["CTRLSHUF_giant"], SVR + " best_protein_q.CTRLSHUF_giant")
put("raw_rev_psms", r6["controls_giant_copies"]["reversed giant"]["psms_q01"], SR + " controls_giant_copies['reversed giant'].psms_q01")
put("raw_shuf_psms", r6["controls_giant_copies"]["shuffled giant"]["psms_q01"], SR + " controls_giant_copies['shuffled giant'].psms_q01")
for key, lab in (("eftu", "elongation factor Tu"), ("l7", "ribosomal protein L7/L12")):
    v = r6["normal_proteins"][lab]
    put(f"raw_{key}_peptides", v["unique_peptides_q01"], SR + f" normal_proteins['{lab}'].unique_peptides_q01")
    put(f"raw_{key}_cov_pct", round(100 * v["coverage_by_unique_frac"]), SR + f" normal_proteins['{lab}'].coverage_by_unique_frac x 100")
put("raw_both", vr["peptides_in_both_searches"], SVR + " peptides_in_both_searches")
put("raw_only_raw", vr["peptides_only_raw"], SVR + " peptides_only_raw")
put("raw_only_mgf", vr["peptides_only_mgf"], SVR + " peptides_only_mgf")
put("raw_mgf_only_after6000", vr["mgf_only_after_6000"], SVR + " mgf_only_after_6000 (start)")
put("raw_in_or_after_repeats", vr["raw_peptides_inside_or_after_repeats"], SVR + " raw_peptides_inside_or_after_repeats")
put("raw_wall_min", 56, "proteomics/search/sage_raw.log 'finished in 3337s' (rounded to min)")
put("prot_peak_last_q", 0.009635, "proteomics/results/full6_giant_peptides.tsv peptide_q of LNTIVTQNK (start 16,949)")

# Sycon ciliatum PXD060105 re-search
SY = os.path.join(ROOT, "proteomics", "sycon")
sy = json.load(open(os.path.join(SY, "results", "summary.json")))
SS = "proteomics/sycon/results/summary.json"
put("sy_ms2", sy["spectra_searched_ms2"], SS + " spectra_searched_ms2")
put("sy_raw_files", sy["raw_files"], SS + " raw_files")
put("sy_db_kept", sy["database_entries_searched"] - 3, SS + " database_entries_searched minus the model and two controls")
put("sy_db_authors", sy["database_entries_authors"], SS + " database_entries_authors (authors' database)")
drop = list(csv.DictReader(open(os.path.join(SY, "db", "dropped_fragments.tsv")), delimiter="\t"))
put("sy_fragments", len(drop), "proteomics/sycon/db/dropped_fragments.tsv rows (fragments of XP_065194412.1)")
put("sy_fragment_max_len", max(int(r["length"]) for r in drop), "same, longest fragment (CAH0927406.1)")
put("sy_psms_spec", sy["spectrum_level"]["target_psms"], SS + " spectrum_level.target_psms")
put("sy_peptides_spec", sy["spectrum_level"]["target_peptides"], SS + " spectrum_level.target_peptides")
put("sy_min_pep_q", round(sy["min_peptide_q_any_target"], 3), SS + " min_peptide_q_any_target")
put("sy_min_prot_q", round(sy["min_protein_q_any_target"], 3), SS + " min_protein_q_any_target")
put("sy_authors_n", sy["authors_reported_proteins"], SS + " authors_reported_proteins")
put("sy_authors_found", sy["spectrum_level"]["authors_proteins_with_psm"], SS + " spectrum_level.authors_proteins_with_psm")
put("sy_giant_psms_spec", sy["spectrum_level"]["giant"]["psms_spectrum_q01"], SS + " spectrum_level.giant.psms_spectrum_q01")
put("sy_giant_psms_any", sy["spectrum_level"]["giant"]["psms_any_q"], SS + " spectrum_level.giant.psms_any_q")
put("sy_rev_psms_any", sy["spectrum_level"]["CTRLREV_giant"]["psms_any_q"], SS + " spectrum_level.CTRLREV_giant.psms_any_q")
put("sy_shuf_psms_any", sy["spectrum_level"]["CTRLSHUF_giant"]["psms_any_q"], SS + " spectrum_level.CTRLSHUF_giant.psms_any_q")
put("sy_calc3_psms", sy["spectrum_level"]["comparison_proteins"]["CAH0794662"]["psms_spectrum_q01"], SS + " spectrum_level.comparison_proteins.CAH0794662")
put("sy_tub_psms", sy["spectrum_level"]["comparison_proteins"]["CAH0788721"]["psms_spectrum_q01"], SS + " spectrum_level.comparison_proteins.CAH0788721")
put("sy_ef2_psms", sy["spectrum_level"]["comparison_proteins"]["CAH0809064"]["psms_spectrum_q01"], SS + " spectrum_level.comparison_proteins.CAH0809064")

am = os.path.join(PR, "authors_mzid_check.json")
if os.path.exists(am):
    a = json.load(open(am))
    put("prot_authors_file", a["file"], "proteomics/results/authors_mzid_check.json file")
    put("prot_authors_len", a["giant_length_in_authors_db"], "proteomics/results/authors_mzid_check.json giant_length_in_authors_db")
    put("prot_authors_peptides", a["distinct_peptides"], "proteomics/results/authors_mzid_check.json distinct_peptides")
    put("prot_authors_fifths", a["starts_per_fifth"], "proteomics/results/authors_mzid_check.json starts_per_fifth")
    put("prot_authors_last_start", a["last_start"], "proteomics/results/authors_mzid_check.json last_start")
    put("prot_authors_pass_peptides", a["passThreshold_only"]["distinct_peptides"], "proteomics/results/authors_mzid_check.json passThreshold_only.distinct_peptides")
    put("prot_authors_pass_fifths", a["passThreshold_only"]["starts_per_fifth"], "proteomics/results/authors_mzid_check.json passThreshold_only.starts_per_fifth")
av = {r["organism"]: r["pride_projects_name_match"] for r in tsv("species_data_availability.tsv")}
for key, org in (("omni", "Omnitrophica bacterium"), ("noc881", "Nocardia sp. NBC_00881"), ("noc1730", "Nocardia sp. NBC_01730"),
                 ("sycon", "Sycon ciliatum"), ("velam", "Candidatus Velamenicoccus archaeovorus")):
    put(f"pride_{key}", int(av[org]), f"data/species_data_availability.tsv pride_projects_name_match, {org}")
pfl = json.load(open(os.path.join(PR, "pride_files_PXD025008.json")))
raw = [x for x in pfl if (x.get("fileCategory") or {}).get("value") == "RAW"]
put("prot_raw_files", len(raw), "proteomics/results/pride_files_PXD025008.json fileCategory RAW (PRIDE API, 13 Nov 2026)")
put("prot_raw_gb", round(sum(x["fileSizeBytes"] for x in raw) / 1e9, 1), "same, sum of fileSizeBytes / 1e9")
put("prot_closed_genome_proteins", sum(1 for l in open(os.path.join(ROOT, "proteomics", "db", "CP019384.1_cds_aa.fasta")) if l.startswith(">")),
    "proteomics/db/CP019384.1_cds_aa.fasta entries (GenBank CDS translations of the closed genome)")
put("raw_giant_last_pct_round", hu(vr["raw_last_residue_pct"]), SVR + " raw_last_residue_pct, rounded")
put("rep_untestable_noruns", sum(1 for r in main if r["md5"] in REP_LOCI and REP_LOCI[r["md5"]] not in rep_job_sp),
    "repeat-flag models in bundle_rep_rna/loci whose species key has no runs (no RNA-seq in ENA)")
put("rep_untestable_nolocus", sum(1 for r in main if "repeat flag" in r["grade_reason_stage1"] and r["md5"] not in REP_LOCI
                                  and r["rnaseq_verdict"].startswith("not-tested (locus not resolved")),
    "repeat-flag models without a locus (not in bundle_rep_rna/loci)")
tw = next(r for r in sens if r["variant"].startswith("record twins"))
put("twin_fix_n", int(tw["sequences_changing_grade"]), "controls/sensitivity/sensitivity_stage12.tsv 'record twins' sequences_changing_grade")
put("twin_fix_big", int(tw["of_which_above_45212"]), "same, of_which_above_45212")
by_org = Counter(r["organism"] for r in main if r["md5"] in REP_LOCI and REP_LOCI[r["md5"]] not in rep_job_sp)
put("rep_noruns_plectranthus", by_org["Plectranthus ornatus"], "repeat-flag models without runs, Plectranthus ornatus")
put("rep_noruns_merlangius", by_org["Merlangius merlangus"], "same, Merlangius merlangus")
put("rep_noruns_coilia", by_org["Coilia grayii"], "same, Coilia grayii")
put("raw_giant_unc_pct", hu(100 - vr["raw_last_residue_pct"]), "100 - raw_last_residue_pct (share of the length after the last peptide)")
lcg = next(r for r in main if r["rep_accession"] == "WP_128698900.1")
put("prot_giant_grade", lcg["grade"], "data/locus_checks.tsv WP_128698900.1 grade")

lb = tsv("leaderboard.tsv")
graded = [r for r in lb if r["grade"] != "X"]
put("lb_loci", len(graded), "data/leaderboard.tsv rows with grade A-D")
put("lb_excluded_rows", len(lb) - len(graded), "data/leaderboard.tsv rows with grade X")
c = Counter(r["grade"] for r in graded)
for gr in "ABCD":
    put(f"lb_{gr}", c[gr], "data/leaderboard.tsv grade, loci")
cb = Counter(r["grade"] for r in graded if int(r["length_aa"]) > PKZ)
put("lb_big", sum(cb.values()), "data/leaderboard.tsv loci with length_aa > 45,212")
for gr in "ABCD":
    put(f"lb_big_{gr}", cb[gr], "data/leaderboard.tsv grade, loci > 45,212 aa")
merged = sum(len([x for x in r["same_gene_other_records"].split("; ") if x]) for r in lb)
put("lb_merged_records", merged, "data/leaderboard.tsv same_gene_other_records, entries")
a_loci = [r for r in graded if r["grade"] == "A"]
put("lb_A_titin_loci", sum(1 for r in a_loci if r["superkingdom"] == "Eukaryota"), "data/leaderboard.tsv grade A, Eukaryota (all titins)")
ea = [int(r["length_aa"]) for r in a_loci if r["superkingdom"] == "Eukaryota"]
put("lb_A_titin_min", min(ea), "same, min length_aa")
put("lb_A_titin_max", max(ea), "same, max length_aa")

# ---- RNA-seq stage -----------------------------------------------------------------------
runs = tsv("rnaseq_runs.tsv")
nosat = {(r["species_key"][:-6], r["run"]) for r in runs if r["species_key"].endswith("_nosat")}
keep = [r for r in runs if (r["species_key"], r["run"]) not in nosat]


def base(s):
    for x in ("_nosat", "_extra"):
        if s.endswith(x):
            return s[:-len(x)]
    return s


put("rna_runs_total", len(runs), "data/rnaseq_runs.tsv rows (all run results)")
put("rna_runs_kept", len(keep), "data/rnaseq_runs.tsv minus runs superseded by a _nosat re-run")
put("rna_runs_superseded", len(runs) - len(keep), "data/rnaseq_runs.tsv runs re-run without the saturation stop")
put("rna_runs_followup", sum(1 for r in runs if r["species_key"].endswith(("_nosat", "_extra"))), "data/rnaseq_runs.tsv species_key ending _nosat/_extra")
put("rna_species", len({base(r["species_key"]) for r in keep}), "data/rnaseq_runs.tsv distinct species (suffixes removed)")
put("rna_runs_v1", sum(1 for r in keep if r["driver"] == "v1"), "data/rnaseq_runs.tsv driver = v1")
v2 = [r for r in keep if r["driver"] == "v2"]
put("rna_runs_v2", len(v2), "data/rnaseq_runs.tsv driver = v2 (kept)")
put("rna_reads_billion", round(sum(int(r["reads_processed"]) for r in v2) / 1e9, 1), "data/rnaseq_runs.tsv sum reads_processed, v2 kept runs")
err = [r for r in keep if "rc=" in r["stop_rule"]]
put("rna_runs_stream_error", len(err), "data/rnaseq_runs.tsv stop_rule contains 'stream error curl rc=' (kept runs)")
put("rna_runs_stream_error_first", sum(1 for r in err if not r["species_key"].endswith(("_nosat", "_extra"))), "same, first RNA box only")
fe = [float(r["frac_of_run_read"]) for r in err if r["frac_of_run_read"]]
put("rna_stream_error_frac_median_pct", round(100 * statistics.median(fe)), "data/rnaseq_runs.tsv frac_of_run_read median over stream-error runs")
put("rna_stream_error_frac_min_pct", round(100 * min(fe)), "same, minimum")
put("rna_runs_saturation", sum(1 for r in keep if "saturation" in r["stop_rule"]), "data/rnaseq_runs.tsv stop_rule contains saturation (kept)")
put("rna_runs_readcap", sum(1 for r in keep if "read cap" in r["stop_rule"]), "data/rnaseq_runs.tsv stop_rule contains read cap (kept)")
fr = [float(r["frac_of_run_read"]) for r in v2 if r["frac_of_run_read"]]
put("rna_frac_read_median_pct", round(100 * statistics.median(fr)), "data/rnaseq_runs.tsv frac_of_run_read median, v2 kept runs")
put("rna_runs_lt50", sum(1 for x in fr if x < 0.5), "same, runs with frac_of_run_read < 0.5")
rs = {r["rep_accession"]: r for r in tsv("rnaseq_support.tsv")}
rv = Counter(r["rnaseq_verdict"] for r in main)
put("rna_full", rv["full"], "data/locus_checks.tsv rnaseq_verdict = full (510)")
put("rna_partial", rv["partial"], "data/locus_checks.tsv rnaseq_verdict = partial")
put("rna_weak", sum(v for k, v in rv.items() if k.startswith("weak")), "data/locus_checks.tsv rnaseq_verdict starts weak")
put("rna_noexpr", rv["no-expression"], "data/locus_checks.tsv rnaseq_verdict = no-expression")
put("rna_noruns", sum(v for k, v in rv.items() if k.startswith("no-rnaseq")), "data/locus_checks.tsv rnaseq_verdict = no-rnaseq-runs")
put("rna_nottested", sum(v for k, v in rv.items() if k.startswith("not-tested")), "data/locus_checks.tsv rnaseq_verdict = not-tested (...)")
put("rna_notselected", rv["not selected"], "data/locus_checks.tsv rnaseq_verdict = not selected (Stage 1 A or D, not a control or fusion check)")
put("rna_prok", rv["not applicable"], "data/locus_checks.tsv rnaseq_verdict = not applicable (prokaryotes)")
put("rna_contradicted", rv["contradicted"], "data/locus_checks.tsv rnaseq_verdict = contradicted")
put("rna_tested_models", sum(1 for r in main if r["rnaseq_verdict"] in ("full", "partial", "no-expression") or r["rnaseq_verdict"].startswith("weak")), "data/locus_checks.tsv models with an RNA verdict from reads")
for key, acc in (("titin", "NP_001254479.2"), ("pkz", "KAL1508129.1"), ("sycon", "XP_065194412.1"), ("halichondria", "XP_064389544.1")):
    r = rs[acc]
    for f in ("runs", "reads_unique", "n_introns", "introns_unique_support", "introns_multi_only", "introns_counted_supported",
              "cds_breadth", "testable_introns", "testable_introns_unsupported"):
        put(f"{key}_{f}", r[f], f"data/rnaseq_support.tsv {acc} column {f}")
put("titin_introns_none", int(rs["NP_001254479.2"]["introns_none"]), "data/rnaseq_support.tsv NP_001254479.2 introns_none")
for key, acc in (("microcaecilia", "XP_030066708.1"), ("rana", "XP_040214433.1"), ("solea", "KAG7502825.1"), ("pantherophis", "XP_087738184.1")):
    r = rs[acc]
    put(f"{key}_testable", r["testable_introns"], f"data/rnaseq_support.tsv {acc} testable_introns")
    put(f"{key}_testable_unsup", r["testable_introns_unsupported"], f"data/rnaseq_support.tsv {acc} testable_introns_unsupported")
    put(f"{key}_verdict", r["verdict"], f"data/rnaseq_support.tsv {acc} verdict")
lift = json.load(open(os.path.join(D, "followups", "relocation_liftover.json")))
put("liftover_moved", lift[next(iter(lift))]["introns_moved"], "data/followups/relocation_liftover.json introns_moved (CAM4536479.1)")
put("liftover_n_introns", len(lift[next(iter(lift))]["introns_lifted"]), "data/followups/relocation_liftover.json introns_lifted length")
r = rs["CAM4536479.1"]
put("cam_after_sup", r["introns_counted_supported"], "data/rnaseq_support.tsv CAM4536479.1 introns_counted_supported (lifted)")
put("cam_after_testable_unsup", r["testable_introns_unsupported"], "data/rnaseq_support.tsv CAM4536479.1 testable_introns_unsupported")

# ---- examples ----------------------------------------------------------------------------
L = {r["rep_accession"]: r for r in lc}
k = L["KAI4551101.1"]
for f in ("length", "n_exons", "genomic_span_bp", "protein_period_aa", "protein_period_identity", "rep25_frac", "gc_cds",
          "dna_tandem_lag_bp", "dna_tandem_identity", "fickett_median", "n_other_frames_open", "min_stops_per_kcodon_other_frames"):
    put(f"ovis_{f}", k[f], f"data/locus_checks.tsv KAI4551101.1 column {f}")
ovr = [r for r in runs if r["species_key"].startswith("Ovis_ammon_polii")]
put("ovis_rna_reads", sum(int(r["reads_processed"]) for r in ovr), "data/rnaseq_runs.tsv Ovis_ammon_polii* reads_processed")
put("ovis_rna_candidates", sum(int(r["candidate_reads"]) for r in ovr), "data/rnaseq_runs.tsv Ovis_ammon_polii* candidate_reads")
ovis_like = sum(1 for r in main if r["organism"].startswith("Ovis ammon polii") and r["protein_period_aa"] == "22")
put("ovis_period22_entries", ovis_like, "data/locus_checks.tsv organism Ovis ammon polii x Ovis aries with protein_period_aa = 22")

for key, acc in (("cirrhinus", "KAK2904762.1"), ("bagarius", "TSR51413.1")):
    x = L[acc]
    put(f"{key}_len", x["length"], f"data/locus_checks.tsv {acc} length")
    put(f"{key}_exons", x["n_exons"], f"data/locus_checks.tsv {acc} n_exons")
    r = rs[acc]
    for f in ("runs", "reads_unique", "n_introns", "introns_counted_supported"):
        put(f"{key}_{f}", r[f], f"data/rnaseq_support.tsv {acc} {f}")
fp = [r for r in tsv("followups/fish_titin_pair.tsv")]
for key, acc in (("cirrhinus", "KAK2904762.1"), ("bagarius", "TSR51413.1")):
    w = [r for r in fp if r["model"] == acc]
    # the second titin copy starts where the matched zebrafish position jumps back to the N-terminus
    pos = [(int(r["window_start_aa"]), r["best_zebrafish_paralogue"],
            int((r["ttna_pos"] if r["best_zebrafish_paralogue"] == "ttna" else r["ttnb_pos"]).split("-")[0])) for r in w]
    jump = next(p for i, p in enumerate(pos) if i and p[2] < pos[i - 1][2] - 10000)
    put(f"{key}_second_copy_start_aa", jump[0], f"data/followups/fish_titin_pair.tsv {acc}: first window whose zebrafish match restarts near the N-terminus")
    first = [p for p in pos if p[0] < jump[0]]
    second = [p for p in pos if p[0] >= jump[0]]
    put(f"{key}_first_ttna", f"{sum(1 for p in first if p[1] == 'ttna')} of {len(first)}", f"data/followups/fish_titin_pair.tsv {acc} windows before the restart with ttna best")
    put(f"{key}_second_ttnb", f"{sum(1 for p in second if p[1] == 'ttnb')} of {len(second)}", f"data/followups/fish_titin_pair.tsv {acc} windows from the restart with ttnb best")

# Cirrhinus per-intron junction reads across the join (RNA results archive)
def intron_table(species_prefix, acc):
    tot, exons = {}, None
    with tarfile.open(os.path.join(D, "rnaseq_run_results.tar.gz")) as t:
        for m in t.getmembers():
            if m.name.startswith("results/" + species_prefix) and m.name.endswith(".json") \
                    and not os.path.basename(m.name).startswith("._"):
                d = json.load(t.extractfile(m))
                for Lx in d.get("loci", []):
                    if Lx["rep_accession"] == acc:
                        exons = exons or Lx["exons"]
                        for i in Lx["introns"]:
                            tot[(i[0], i[1])] = tot.get((i[0], i[1]), 0) + i[2]
    return tot, exons


tot, exons = intron_table("Cirrhinus_molitorella", "KAK2904762.1")
strand = L["KAK2904762.1"]["strand"]
ex = sorted((e[0], e[1]) for e in exons)
if strand == "-":
    ex = ex[::-1]
cum, aa_at = 0, []
order = []
for i, (a, b) in enumerate(ex[:-1]):
    cum += b - a + 1
    nxt = ex[i + 1]
    key = (b + 1, nxt[0] - 1) if strand == "+" else (nxt[1] + 1, a - 1)
    order.append((cum // 3, tot.get(key, 0)))
jidx = min(range(len(order)), key=lambda i: abs(order[i][0] - 30675))
put("cirrhinus_join_intron_aa", order[jidx][0], "rnaseq_run_results.tar.gz Cirrhinus results: intron nearest aa 30,675 (CDS position / 3)")
put("cirrhinus_join_intron_reads", order[jidx][1], "same: unique exact junction reads summed over runs")
put("cirrhinus_join_before_median", statistics.median([x[1] for x in order[jidx - 5:jidx]]), "same: median of the 5 introns before")
put("cirrhinus_join_after_median", statistics.median([x[1] for x in order[jidx + 1:jidx + 6]]), "same: median of the 5 introns after")

put("ovis_entries", sum(1 for r in main if r["organism"].startswith("Ovis ammon polii")), "data/locus_checks.tsv organism Ovis ammon polii x Ovis aries")
put("ovis_period_multiple22", sum(1 for r in main if r["organism"].startswith("Ovis ammon polii") and r["protein_period_aa"]
                                 and int(r["protein_period_aa"]) % 22 == 0), "same, protein_period_aa a multiple of 22")
syc = next(r for r in tsv("ncbi_records_parsed.tsv") if r["accession"] == "XP_065194412.1")
put("sycon_gap_patch_bp", int(re.search(r"added (\d+) bases", syc["model_modifications"]).group(1)),
    "data/ncbi_records_parsed.tsv XP_065194412.1 model_modifications")
put("sycon_best_homolog_len", L["XP_065194412.1"]["best_other_species_len"], "data/locus_checks.tsv XP_065194412.1 best_other_species_len")
put("sycon_len", L["XP_065194412.1"]["length"], "data/locus_checks.tsv XP_065194412.1 length")
for key, acc in (("nbc00881", "WP_406270362.1"), ("nbc01730", "WP_327112628.1")):
    put(f"{key}_len", L[acc]["length"], f"data/locus_checks.tsv {acc} length")
    put(f"{key}_ratio", L[acc]["length_ratio_query_over_best"], f"data/locus_checks.tsv {acc} length_ratio_query_over_best")
    put(f"{key}_best_len", L[acc]["best_other_species_len"], f"data/locus_checks.tsv {acc} best_other_species_len")
    put(f"{key}_orf_bp", L[acc]["cds_bp"], f"data/locus_checks.tsv {acc} cds_bp")
put("pkz_uniprot_pe", next(r["Protein existence"] for r in up if r["Entry"] == "A0AB34IYJ6"),
    "data/uniprotkb_ge35992.tsv A0AB34IYJ6 'Protein existence'")
put("stage1_cpu_hours", 4.4, "notes/FEASIBILITY.md section 11 (measured CPU time of Stage 1; not in a data file)")

# stops per 1,000 codons in the six frames of the KAI4551101.1 CDS (full join from GenPept)
gp = gzip.open(os.path.join(D, "raw", "ncbi_gp.txt.gz"), "rt").read()
i0 = gp.find("/coded_by=", gp.find("KAI4551101"))
coded_by = re.sub(r"\s+", "", gp[i0 + 11:gp.find('"', i0 + 11)])
lseq = "".join(l.strip().upper() for l in gzip.open(os.path.join(D, "manual", "KAI4551101_locus.fa.gz"), "rt") if l[0] != ">")
hdr = gzip.open(os.path.join(D, "manual", "KAI4551101_locus.fa.gz"), "rt").readline()
off = int(hdr.split(":")[1].split("-")[0])
parts = [(int(a), int(b)) for a, b in re.findall(r":<?(\d+)\.\.>?(\d+)", coded_by)]
cds = "".join(lseq[a - off:b - off + 1] for a, b in sorted(parts))
COMP = str.maketrans("ACGTN", "TGCAN")
if coded_by.startswith("complement"):
    cds = cds.translate(COMP)[::-1]
fr = []
for strand in (cds, cds.translate(COMP)[::-1]):
    for f0 in range(3):
        cod = [strand[i:i + 3] for i in range(f0, len(strand) - 2, 3)]
        fr.append(round(1000 * sum(c in ("TAA", "TAG", "TGA") for c in cod) / len(cod), 1))
put("ovis_stops_six_frames", fr, "computed from data/manual/KAI4551101_locus.fa.gz + CDS join in data/raw/ncbi_gp.txt.gz (+1 = annotated)")

# ---- DNA checks --------------------------------------------------------------------------
dna = {}
for f in os.listdir(os.path.join(D, "dna_checks")):
    j = json.load(open(os.path.join(D, "dna_checks", f)))
    dna[j["job_id"]] = j
SHORT = {"omnitrophota_DAOOUB_HOW35280": "oub", "omnitrophota_DAOPMG_HPD84636": "pmg", "omnitrophota_DAOTDT_HRZ03494": "tdt",
         "nocardia_NBC_01730_WSG58785": "nbc01730", "nocardia_NBC_00881_WSY60125": "nbc00881"}
for tag, j in dna.items():
    il = j["illumina"]
    name = SHORT[tag]
    src = f"data/dna_checks/{j['job_id']}.json"
    put(f"dna_{name}_reads", j["illumina_reads_processed"], src + " illumina_reads_processed")
    put(f"dna_{name}_reads_expected", j["illumina_expected_reads"], src + " illumina_expected_reads")
    put(f"dna_{name}_depth", il["gene_median_depth_unique_windows"], src + " illumina.gene_median_depth_unique_windows")
    put(f"dna_{name}_ratio", il["depth_ratio_gene_vs_contig"], src + " illumina.depth_ratio_gene_vs_contig")
    put(f"dna_{name}_span_median", il["median_spanning_per_boundary"], src + " illumina.median_spanning_per_boundary")
    put(f"dna_{name}_span_min_testable", il["min_spanning_unique_boundaries"], src + " illumina.min_spanning_unique_boundaries")
    put(f"dna_{name}_discordant_pct", round(100 * il["max_discordant_frac_window"], 1), src + " illumina.max_discordant_frac_window")
    put(f"dna_{name}_insert", il["median_insert"], src + " illumina.median_insert")
    put(f"dna_{name}_untestable", len(il["boundaries_untestable_for_library"]), src + " illumina.boundaries_untestable_for_library")
    put(f"dna_{name}_boundaries", il["boundaries"], src + " illumina.boundaries")
    put(f"dna_{name}_verdict", j["final_verdict"], src + " final_verdict")
    put(f"dna_{name}_control", j["negative_control_illumina"]["control_detected_the_join"], src + " negative_control_illumina.control_detected_the_join")
    put(f"dna_{name}_v2", (j.get("previous") or {}).get("final_verdict", "n/a (judged under v3 only)"), src + " previous.final_verdict")
    E = il.get("expected_spanning_per_boundary") or []
    S = il["spanning_per_boundary"]
    unt = [k for k in range(len(S)) if E and E[k] is not None and E[k] < 5]
    put(f"dna_{name}_E_at_untestable", [E[k] for k in unt], src + " expected_spanning_per_boundary at untestable boundaries")
    put(f"dna_{name}_obs_at_untestable", [S[k] for k in unt], src + " spanning_per_boundary at untestable boundaries")
    tst = sorted((S[k], E[k]) for k in range(len(S)) if E and E[k] is not None and E[k] >= 5)
    if tst:
        put(f"dna_{name}_weakest", tst[0], src + " weakest testable boundary (observed, expected)")
    if "ont" in j:
        put(f"dna_{name}_ont_reads", j["ont_reads_processed"], src + " ont_reads_processed")
        sb = j["ont"].get("short_read_repeat_blocks_spanned_by_long_reads", [])
        put(f"dna_{name}_ont_span_blocks", [b["ont_reads_spanning"] for b in sb], src + " ont.short_read_repeat_blocks_spanned_by_long_reads")

# internal repeat of the 85,803-aa ORF, from the committed MAG (data/manual/DAOTDT01.fasta.gz)
seqs, name, buf = {}, None, []
with gzip.open(os.path.join(D, "manual", "DAOTDT01.fasta.gz"), "rt") as f:
    for line in f:
        if line[0] == ">":
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0].split("|")[-1], []
        else:
            buf.append(line.strip().upper())
seqs[name] = "".join(buf)
tdt = dna["omnitrophota_DAOTDT_HRZ03494"]
gs, ge = tdt["gene"]
gene = seqs[tdt["contig"]][gs - 1:ge]
put("orf_bp", len(gene), "data/manual/DAOTDT01.fasta.gz, CDS of HRZ03494 (contig DAOTDT010000001.1, gene coordinates from the DNA-check JSON)")
put("contig_len_tdt", len(seqs[tdt["contig"]]), "data/manual/DAOTDT01.fasta.gz contig DAOTDT010000001.1 length")
K = 40
idx = {}
for i in range(len(gene) - K + 1):
    idx.setdefault(gene[i:i + K], []).append(i)
pairs = Counter()
for v in idx.values():
    if 1 < len(v) <= 3:
        for a in v:
            for b in v:
                if b > a + 1000:
                    pairs[b - a] += 1
lag, _ = pairs.most_common(1)[0]
pos = sorted({a for v in idx.values() if 1 < len(v) <= 3 for a in v for b in v if b - a == lag})
a0 = pos[0]
a1 = pos[-1] + K
copy_a, copy_b = gene[a0:a1], gene[a0 + lag:a1 + lag]
put("rep_lag_bp", lag, "computed: most common offset of repeated 40-mers in the CDS (paper/make_numbers.py)")
put("rep_lag_codons", lag / 3, "rep_lag_bp / 3")
put("rep_copy_bp", a1 - a0, "computed: span of 40-mers shared at that offset")
put("rep_copy_a", (a0 + 1, a1), "computed: copy A, CDS bp (1-based from the A of ATG)")
put("rep_copy_b", (a0 + lag + 1, a1 + lag), "computed: copy B, CDS bp")
put("rep_diffs", sum(1 for x, y in zip(copy_a, copy_b) if x != y), "computed: mismatches between the two copies (ungapped)")
put("orf_len_aa", (N["orf_bp"][0] - 3) // 3, "orf_bp minus the stop codon, / 3")
put("rep_copy_bp_kb", round(N["rep_copy_bp"][0] / 1000, 2), "rep_copy_bp / 1000")
put("dna_oub_weakest_obs", N["dna_oub_weakest"][0][0], "dna_oub_weakest observed pairs")
put("dna_oub_weakest_exp", hu(N["dna_oub_weakest"][0][1]), "dna_oub_weakest expected pairs, rounded")
put("orf_len_if_excised", (len(gene) - lag - 3) // 3, "(orf_bp - rep_lag_bp - stop codon) / 3")

# the three copies of the 85,803-aa ORF: protein identity (unique sequences file) and contigs
prot, name, buf = {}, None, []
with gzip.open(os.path.join(D, "fasta", "unique_ge35992.fasta.gz"), "rt") as f:
    for line in f:
        if line[0] == ">":
            if name:
                prot[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
prot[name] = "".join(buf)
md5 = {r["rep_accession"]: r["md5"] for r in lc}
cp = ["UPI002B62E083"] + L["UPI002B62E083"]["near_identical_copy_accessions"].split(";")
ids = []
for i in range(len(cp)):
    for j in range(i + 1, len(cp)):
        a, b = prot[md5[cp[i]]], prot[md5[cp[j]]]
        ids.append(sum(x == y for x, y in zip(a, b)) / max(len(a), len(b)))
put("omni_copies", len(cp), "data/locus_checks.tsv UPI002B62E083 near_identical_copy_accessions + itself")
put("omni_aa_identity_min_pct", round(100 * min(ids), 2), "data/fasta/unique_ge35992.fasta.gz, ungapped identity of the 3 copies (min)")
put("omni_aa_identity_max_pct", round(100 * max(ids), 2), "same (max)")
put("omni_contig_lengths", [dna[k]["contig_length"] for k in ("omnitrophota_DAOTDT_HRZ03494", "omnitrophota_DAOPMG_HPD84636",
                                                              "omnitrophota_DAOOUB_HOW35280")], "data/dna_checks/*.json contig_length")
put("omni_n_in_cds", L["UPI002B62E083"]["N_in_cds"], "data/locus_checks.tsv UPI002B62E083 N_in_cds")
put("hmm_panel", sum(1 for l in open(os.path.join(D, "hmm", "panel.hmm")) if l.startswith("NAME")), "data/hmm/panel.hmm NAME lines")
a2 = [r for r in main if r["grade_stage2_rna"] == "A"]
put("s2_A_titin", sum(1 for r in a2 if r["best_other_species_name"] == "Titin"), "data/locus_checks.tsv grade_stage2_rna A with best_other_species_name Titin")
put("s1_fix_B_to_C", 48, "notes/FEASIBILITY.md section 14 and notes/RESULTS.md section 2 (Stage 1 B 267 -> 219, C 53 -> 101 after the S2 fix)")
put("titin_downgrades_withdrawn", 5, "notes/RESULTS.md section 3 ('moved five titin models from B to C', withdrawn)")

put("cirrhinus_first_end", int(N["cirrhinus_second_copy_start_aa"][0]) - 1, "cirrhinus_second_copy_start_aa - 1")
json.dump({k: v for k, (v, s) in N.items()}, open(os.path.join(ROOT, "paper", "numbers.json"), "w"), indent=1, default=str)
with open(os.path.join(ROOT, "paper", "NUMBERS_generated.tsv"), "w") as f:
    f.write("key\tvalue\tsource\n")
    for k, (v, s) in N.items():
        f.write(f"{k}\t{v}\t{s}\n")


def tex_value(v):
    """How a value is printed in the paper: thousands separators, q-values to 1-2 significant digits,
    a pair as a range with an en dash, a list joined by commas."""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, int):
        return f"{v:,}"
    if isinstance(v, float):
        if v == int(v) and abs(v) >= 1:
            return f"{int(v):,}"
        if abs(v) < 0.01:
            return f"{v:.2g}"
        if abs(v) < 1 and len(repr(v)) > 6:
            return f"{v:.2g}"
        return f"{v:,}"
    if isinstance(v, (list, tuple)):
        if len(v) == 2 and all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in v):
            return tex_value(v[0]) + "--" + tex_value(v[1])
        return ", ".join(tex_value(x[0]) if isinstance(x, (list, tuple)) else tex_value(x) for x in v)
    t = str(v)
    if re.fullmatch(r"\d{4,}", t):                 # integers read from TSV files as text
        return f"{int(t):,}"
    for a, b in (("\\", "\\textbackslash{}"), ("&", "\\&"), ("%", "\\%"), ("#", "\\#"), ("_", "\\_"),
                 ("$", "\\$"), (">=", "$\\geq$"), ("<=", "$\\leq$")):
        t = t.replace(a, b)
    return t


with open(os.path.join(ROOT, "paper", "numbers.tex"), "w") as f:
    f.write("% Generated by paper/make_numbers.py. Do not edit. Sources: paper/NUMBERS_generated.tsv\n")
    for k, (v, s) in N.items():
        if k.startswith("_"):
            continue
        f.write(f"\\nmdef{{{k}}}{{{tex_value(v)}}}\n")
for k, (v, s) in N.items():
    print(f"{k:34s} {v}")
