"""Job 1, step 4: merge Stage 0 flags, locus DNA checks (11), orthology checks (12) and
RefSeq model evidence into data/locus_checks.tsv with a provisional grade (Stage 1,
before any new RNA-seq) and the reason in words. Controls: PKZILLA-1 (in the set) and
human titin NP_001254479.2 (added by 10b_controls.py).

Rules (thresholds fixed before looking at the controls' grades; see notes/FEASIBILITY.md):
  Contradictions -> D
    D1 repeat translation: protein self-identity at best period >= 0.6 or rep25 >= 0.8,
       and recognised-domain coverage < 5 %; since 13 Nov (D1 rule v2) a contradiction only
       together with D2 or D3, otherwise a 'repeat flag' (positive controls: controls/d1)
    D2 non-coding DNA: >= 50 % of 300-bp CDS windows have Fickett TESTCODE <= 0.74,
       and domain coverage < 5 %
    D3 frame-independent ORF: >= 2 of the 5 alternative frames also open
       (< 5 stops per 1,000 codons), and domain coverage < 5 %
    D4 assembly: >= 1 % N in the CDS, or the model joins >= 2 contigs
    D5 fusion: >= 2 disjoint segments each covering >= 70 % of a different complete
       reference protein (>= max(1,000 aa, 20 % of the query) and <= 60 % of the query),
       together spanning >= 60 % of the query, and no other-species homologue within
       +-20 % of the query length (not applied to PKS/NRPS, whose
       modules are expected to match smaller complete synthases; flagged instead)
  Structure OK = no D trigger, no N in CDS, gene >= 1 kb from both contig ends, coding
       signal (Fickett median >= 0.95 and <= 20 % non-coding windows, or GC2 lowest in
       >= 80 % of 3-kb windows)
  Independent support lines
    S1 orthologue length: >= 1 other-species UniProtKB protein hit on >= 50 % of the
       query windows with length within +-20 %
    S2 other assembly: near-identical copy (>= 95 % id, >= 90 % windows, length +-5 %)
       of the same species from ANOTHER assembly (contig accession prefix differs;
       isoforms of one gene and same-assembly UniProtKB entries do not count)
    S3 transcripts: NCBI Gnomon model evidence reports >= 1 RNA-seq sample supporting
       all annotated introns
    S4 curated RefSeq record (NP_/NM_)
    S6 peptides (13 Nov): >= 3 unique peptides at 1 % FDR and >= 1 in each of the 5 length
       fifths (peptides along the whole model); peptides in fewer fifths are noted in the
       grade reason but are not a support line (GP_S6_FIFTHS sets the number of fifths)
  A = structure OK + (S3 or S4) + (S1 or S2)   [S6 counts as a support line for B only]
  B = structure OK + any of S1-S4 or S6
  C = no D trigger but not B
"""
import csv
import os

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
# thresholds (defaults as published; environment overrides exist only for the sensitivity
# analysis in controls/sensitivity/run_sensitivity.py)
T_D1_ID = float(os.environ.get("GP_D1_ID", "0.6"))        # D1 protein self-identity at the best period
T_D1_REP25 = float(os.environ.get("GP_D1_REP25", "0.8"))  # D1 share of residues in repeated 25-mers
T_DOM = float(os.environ.get("GP_DOM", "0.05"))           # D1-D3 apply only below this domain coverage
T_D2_NC = float(os.environ.get("GP_D2_NC", "0.5"))        # D2 share of non-coding Fickett windows
T_D3_FRAMES = int(os.environ.get("GP_D3_FRAMES", "2"))    # D3 open alternative frames
T_FICKETT = float(os.environ.get("GP_FICKETT_MED", "0.95"))  # coding signal: Fickett median
# D1 rule version. v1 (used until 13 Nov): a repeat translation alone is a contradiction.
# v2 (13 Nov, after the positive controls in controls/d1: v1 flagged 10 of 176 real repetitive
# proteins with protein-level evidence, among them silks and serine-rich adhesins): the repeat
# trigger counts as a contradiction only together with D2 or D3, i.e. DNA-level evidence that
# the sequence does not behave like a protein-coding gene; alone it is reported as a flag.
D1_RULE = os.environ.get("GP_D1_RULE", "v2")
TWIN_FIX = os.environ.get("GP_TWIN_FIX", "1") == "1"  # 0 = record-twin handling before the 13 Nov fix
T_S6_FIFTHS = int(os.environ.get("GP_S6_FIFTHS", "5"))   # S6 length fifths that need a unique peptide


def load(name, key="md5"):
    p = os.path.join(D, name)
    return {r[key]: r for r in csv.DictReader(open(p), delimiter="\t")} if os.path.exists(p) else {}


s0 = load("stage0_table.tsv")
dna = load("locus_dna_checks.tsv")
orth = load("orthology_checks.tsv")
specs_order = list(dna)              # largest first (specs were written in length order)
CONTROL_TITIN = [m for m in dna if m.startswith("control_")]
ORG_BY_ACC = {r["rep_accession"]: r["organism"] for r in s0.values()}
CONTIG_BY_ACC = {r["rep_accession"]: r.get("contig", "") for r in dna.values()}


def assembly_id(contig):
    """Assembly a contig belongs to: the WGS/INSDC accession prefix (letters), or 'RefSeq'
    for NC_/NW_/NZ_/NT_ sequences (RefSeq annotates one assembly per species)."""
    import re as _re
    if not contig:
        return ""
    if contig.startswith("NZ_") and TWIN_FIX:
        # a RefSeq copy (NZ_) of an INSDC record is the SAME assembly (fixed 13 Nov 2026: three
        # RefSeq/GenBank twin pairs had been counted as 'another assembly' for S2)
        contig = contig[3:]
    elif contig.startswith(("NC_", "NW_", "NT_", "NZ_")):
        return "RefSeq"
    m = _re.match(r"^([A-Z]+)", contig)
    return m.group(1) if m else contig


# S6 (added 13 Nov 2026 with the PXD025008 reanalysis, after its result was seen): unique peptides at
# 1 % FDR (PSM, peptide and protein level) in every fifth of the protein length, i.e. along the whole
# model. A first draft counted >= 3 fifths, which the one tested entry meets with 182 of its 186
# peptides in the first fifth; that reading is kept only as a sensitivity variant (GP_S6_FIFTHS=3).
# data/followups/peptide_evidence.tsv
PEPTIDES = {}
_pp = os.path.join(D, "followups", "peptide_evidence.tsv")
if os.path.exists(_pp):
    PEPTIDES = {r["accession"]: r for r in csv.DictReader(open(_pp), delimiter="\t")}


def species_key(org):
    """Genus + species, ignoring a 'Candidatus' prefix and strain text (fixed 13 Nov 2026: the
    orthology search had counted 'Velamenicoccus archaeovorus', the UniProtKB copy of the same
    genome, as another species for 'Candidatus Velamenicoccus archaeovorus')."""
    import re as _re
    o = _re.sub(r"^(Candidatus|Ca\.)\s+", "", (org or "").replace('"', "").strip())
    o = _re.sub(r"\s*\(.*\)", "", o)
    return " ".join(o.split()[:2]).lower()


def f(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def prefix(acc):
    return acc.split(".")[0][:6] if acc else ""


rows = []
for m in specs_order + [x for x in s0 if x not in dna and s0[x]["triage"] != "T-excl"]:
    a = s0.get(m, {})
    d = dna.get(m, {})
    o = orth.get(m, {})
    ctrl = m.startswith("control_") or a.get("rep_accession") == "KAL1508129.1"
    length = int(a["length"]) if a.get("length") else 35991      # control: NP_001254479.2
    dom = max(f(a.get("panel_cov_frac")), f(a.get("uniparc_pfam_cov_frac")))
    pid = f(a.get("top_period_identity"))
    rep25 = f(a.get("rep25_frac"))
    has_dna = bool(d.get("cds_bp"))
    cds_bp = f(d.get("cds_bp"), 1)
    reasons, contra, support = [], [], []
    d1 = (pid >= T_D1_ID or rep25 >= T_D1_REP25) and dom < T_DOM
    d2 = has_dna and f(d.get("fickett_noncoding_frac")) >= T_D2_NC and dom < T_DOM
    d3 = has_dna and int(d.get("n_other_frames_open") or 0) >= T_D3_FRAMES and dom < T_DOM
    d1_text = f"D1 repeat translation (period {a.get('top_period_aa')} aa at identity {pid:.2f}, rep25 {rep25:.2f}, no domains)"
    if d1 and (D1_RULE == "v1" or d2 or d3):
        contra.append(d1_text)
    elif d1:
        reasons.append("repeat flag (D1 pattern without DNA-level contradiction; real repetitive proteins such as "
                       f"silks show it too): period {a.get('top_period_aa')} aa at identity {pid:.2f}, rep25 {rep25:.2f}")
    if d2:
        contra.append(f"D2 non-coding DNA signal (Fickett: {f(d['fickett_noncoding_frac']):.0%} of windows non-coding)")
    if d3:
        contra.append(f"D3 ORF open in {d['n_other_frames_open']} other frames too")
    if has_dna and (f(d.get("N_in_cds")) >= 0.01 * cds_bp or int(d.get("n_contigs_in_model") or 1) > 1):
        contra.append(f"D4 assembly gap in CDS ({d.get('N_in_cds')} N) or multi-contig model")
    fusion = int(o.get("complete_proteins_on_disjoint_segments") or 0) >= 2 and \
        f(o.get("disjoint_segments_cover_frac")) >= 0.6 and \
        int(o.get("n_other_species_hits_cov50_len_within_20pct") or 0) == 0
    modular = a.get("panel_class", "").startswith(("PKS", "NRPS"))
    if fusion and not modular:
        contra.append("D5 segments match different complete proteins: " + o.get("fusion_segments", "")[:160])
    elif fusion and modular:
        reasons.append("segments match smaller complete synthases (expected for modular PKS/NRPS)")
    coding = has_dna and ((f(d.get("fickett_median")) >= T_FICKETT and f(d.get("fickett_noncoding_frac")) <= 0.2)
                          or f(d.get("coding_windows_gc2_lowest_frac")) >= 0.8)
    near_end = has_dna and d.get("contig_length") not in ("", None) and \
        min(f(d.get("dist_to_contig_start"), 1e9), f(d.get("dist_to_contig_end"), 1e9)) < 1000
    structure_ok = has_dna and not contra and f(d.get("N_in_cds")) == 0 and coding and not near_end
    # support lines
    same_species = TWIN_FIX and species_key(o.get("best_other_species_org", "")) == species_key(a.get("organism", "Homo sapiens"))
    if int(o.get("n_other_species_hits_cov50_len_within_20pct") or 0) >= 1 and not same_species:
        support.append(f"S1 orthologue length consistent ({o.get('best_other_species_org', '')[:30]} "
                       f"{o.get('best_other_species_len')} aa, ratio {o.get('length_ratio_query_over_best')})")
    # S2 counts only copies from ANOTHER ASSEMBLY of the same species: isoforms of one gene
    # (same GeneID / same contig) and same-assembly UniProtKB entries are not independent.
    # (Fixed 12 Nov 2026: the first version also counted isoforms and same-species UniProtKB
    # entries, which can come from the same assembly.)
    my_asm = assembly_id(d.get("contig", "") or a.get("locus_key", ""))
    copies = [c for c in (o.get("near_identical_copy_accessions") or "").split(";")
              if c and ORG_BY_ACC.get(c, "").split(" (")[0] == a.get("organism", "Homo sapiens").split(" (")[0]
              and assembly_id(CONTIG_BY_ACC.get(c, "")) not in ("", my_asm)]
    s2_ev = []
    if copies:
        s2_ev.append("near-identical copy in another assembly: " + ";".join(copies)[:60])
    if s2_ev:
        support.append("S2 " + "; ".join(s2_ev))
    pe = PEPTIDES.get(a.get("rep_accession", ""))
    if pe and int(pe["unique_peptides"]) >= 3 and int(pe["fifths_with_unique_peptides"]) >= T_S6_FIFTHS:
        support.append(f"S6 peptides: {pe['unique_peptides']} unique peptides at 1 % FDR in {pe['fifths_with_unique_peptides']} "
                       f"of 5 length fifths ({pe['dataset']})")
    elif pe:
        reasons.append(f"peptides for part of the model only ({pe['unique_peptides']} unique peptides at 1 % FDR, "
                       f"aa {pe['first_peptide_start']}-{pe['last_peptide_end']}, {pe['fifths_with_unique_peptides']} of 5 length fifths; "
                       f"{pe['dataset']}; not a support line)")
    if int(a.get("rnaseq_samples_all_introns") or 0) >= 1:
        support.append(f"S3 NCBI Gnomon: {a.get('rnaseq_feature_coverage_pct')}% RNA-seq coverage, "
                       f"{a.get('rnaseq_samples_all_introns')} sample(s) support all introns")
    if a.get("rep_accession", "").startswith("NP_") or m.startswith("control_"):
        support.append("S4 curated RefSeq record")
    s_ids = {x[:2] for x in support}
    if contra:
        grade = "D"
    elif structure_ok and ({"S3", "S4"} & s_ids) and ({"S1", "S2"} & s_ids):
        grade = "A"
    elif structure_ok and s_ids:
        grade = "B"
    else:
        grade = "C"
        if not has_dna:
            reasons.append("locus not resolved, no DNA checks")
        elif near_end:
            reasons.append("gene ends < 1 kb from a contig end (possibly truncated or joined at an assembly break)")
        elif not coding:
            reasons.append(f"weak coding signal (Fickett median {d.get('fickett_median')}, "
                           f"non-coding windows {d.get('fickett_noncoding_frac')})")
        if not support:
            reasons.append("no independent support line (no consistent orthologue, no other assembly, no transcript evidence)")
    reason = "; ".join(contra + support + reasons)
    row = {
        "md5": m, "rep_accession": a.get("rep_accession") or d.get("rep_accession"),
        "control": "positive" if ctrl else "",
        "length": length,
        "organism": a.get("organism", "Homo sapiens"), "superkingdom": a.get("superkingdom", "Eukaryota"),
        "databases": a.get("databases", "RefSeq"), "uniprotkb": a.get("uniprotkb", ""), "uniprotkb_PE": a.get("uniprotkb_PE", ""),
        "stage0_triage": a.get("triage", "control"), "panel_class": a.get("panel_class", ""),
        "domain_cov_frac": round(dom, 3),
        "protein_period_aa": a.get("top_period_aa", ""), "protein_period_identity": a.get("top_period_identity", ""),
        "rep25_frac": a.get("rep25_frac", ""), "lowcomp_frac": a.get("lowcomp_frac", ""), "x_frac": a.get("x_frac", ""),
    }
    for c in ("locus_source", "contig", "strand", "n_contigs_in_model", "n_exons", "cds_bp", "mean_exon_bp",
              "median_exon_bp", "max_exon_bp", "mean_intron_bp", "max_intron_bp", "genomic_span_bp",
              "contig_length", "dist_to_contig_start", "dist_to_contig_end", "N_in_cds", "gap_runs_in_span",
              "N_in_span", "translation_identity", "translation_len_diff", "gc_cds", "gc1", "gc2", "gc3",
              "coding_windows_gc2_lowest_frac", "fickett_median", "fickett_coding_frac", "fickett_noncoding_frac",
              "internal_stops_annot_frame", "min_stops_per_kcodon_other_frames",
              "mean_stops_per_kcodon_other_frames", "n_other_frames_open", "dna_tandem_lag_bp",
              "dna_tandem_identity", "locus_note"):
        row[c] = d.get(c, "")
    row["single_orf_or_exon_chain"] = ("" if not d.get("n_exons") else
                                       "single ORF" if d["n_exons"] == "1" else f"exon chain ({d['n_exons']})")
    for c in ("windows", "windows_with_ref_hit", "best_other_species_hit", "best_other_species_org",
              "best_other_species_name", "best_other_species_len", "best_other_species_windows",
              "best_other_species_pid", "length_ratio_query_over_best", "n_other_species_hits_cov50",
              "n_other_species_hits_cov50_len_within_20pct", "best_same_species_other_entry",
              "best_same_species_len", "best_same_species_windows", "identical_entries_in_uniprotkb",
              "identical_entry_accessions",
              "targets_to_cover_80pct_hit_windows", "complete_proteins_on_disjoint_segments",
              "disjoint_segments_cover_frac", "fusion_segments",
              "near_identical_copies_in_giant_set", "near_identical_copy_accessions"):
        row[c] = o.get(c, "")
    row["refseq_rnaseq_coverage_pct"] = a.get("rnaseq_feature_coverage_pct", "")
    row["refseq_samples_all_introns"] = a.get("rnaseq_samples_all_introns", "")
    row["stage0_flags"] = a.get("flags", "")
    row["grade"] = grade
    row["grade_reason"] = reason
    rows.append(row)

rows.sort(key=lambda r: -int(r["length"]))
cols = list(rows[0].keys())
with open(os.path.join(D, "locus_checks.tsv"), "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in rows:
        fo.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
from collections import Counter  # noqa: E402
print("rows", len(rows), Counter(r["grade"] for r in rows))
print(">45212", Counter(r["grade"] for r in rows if int(r["length"]) > 45212))
for r in rows:
    if r["control"]:
        print("control", r["rep_accession"], r["grade"], "|", r["grade_reason"][:200])
