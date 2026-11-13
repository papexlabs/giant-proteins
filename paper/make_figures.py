"""Figures for paper/manuscript.md, made only from files in this repository.

Run: .venv/bin/python paper/make_figures.py   (needs matplotlib + numpy in the project venv)
Writes paper/figures/fig{1..5}.{pdf,png}. Colours: Okabe-Ito (colour-blind safe); grades are
also labelled in text, so no figure depends on colour alone.
"""
import csv
import gzip
import json
import os
import re
import tarfile
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
D = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "paper", "figures")
os.makedirs(OUT, exist_ok=True)
TITIN, PKZ = 35991, 45212
GC = {"A": "#009E73", "B": "#56B4E9", "C": "#E69F00", "D": "#D55E00"}     # Okabe-Ito
BLUE, ORANGE, GREY, BLACK, PURPLE = "#0072B2", "#E69F00", "#999999", "#000000", "#CC79A7"
GRADE_NAME = {"A": "A  supported full-length model", "B": "B  supported in parts", "C": "C  no contradiction, no support line",
              "D": "D  contradicted by a check"}
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 8,
    "axes.titlesize": 8.5, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none",
    "figure.dpi": 150, "savefig.dpi": 300})
COMP = str.maketrans("ACGTN", "TGCAN")


def tsv(name):
    return list(csv.DictReader(open(os.path.join(D, name)), delimiter="\t"))


def save(fig, name):
    for ext in ("pdf", "png"):
        md = {"CreationDate": None} if ext == "pdf" else {}       # byte-stable PDFs across builds
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight", metadata=md)
    plt.close(fig)


def panel(ax, letter):
    ax.text(-0.13, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def fasta_gz(path, key=lambda h: h.split()[0]):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] == ">":
                if name:
                    seqs[name] = "".join(buf)
                name, buf = key(line[1:]), []
            else:
                buf.append(line.strip().upper())
    seqs[name] = "".join(buf)
    return seqs


lc = tsv("locus_checks.tsv")
main = [r for r in lc if not r["md5"].startswith("control_")]

# ---- Figure 1: lengths and grades ----------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 4.6), sharex=True, gridspec_kw={"height_ratios": [1, 1.3]})
bins = np.arange(35000, 101000, 1000)
allL = [int(r["length"]) for r in main]
h, _ = np.histogram(allL, bins=bins)
ax1.bar(bins[:-1], h, width=1000, align="edge", color=GREY, edgecolor="white", linewidth=0.3)
ax1.set_yscale("log")
ax1.set_ylim(0.7, 2000)
ax1.set_yticks([1, 10, 100])
ax1.set_yticklabels(["1", "10", "100"])
ax1.set_ylabel("Sequences per 1,000 aa")
ax1.set_title(f"{len(main)} retained sequences longer than titin (bars: count, log scale; grades in b)", loc="left")
panel(ax1, "a")
rng = np.random.default_rng(1)
rows = {"A": 3, "B": 2, "C": 1, "D": 0}
for gr in "ABCD":
    L = np.array([int(r["length"]) for r in main if r["grade"] == gr])
    y = rows[gr] + rng.uniform(-0.28, 0.28, len(L))
    ax2.scatter(L, y, s=7, color=GC[gr], edgecolor="none", alpha=0.85)
    nb = int((L > PKZ).sum())
    ax2.text(100800, rows[gr], f"{len(L)} ({nb} > PKZILLA-1)", va="center", fontsize=6.5)
ax2.set_yticks([rows[g] for g in "ABCD"])
ax2.set_yticklabels(["A", "B", "C", "D"])
ax2.set_ylabel("Grade")
ax2.set_ylim(-0.6, 3.6)
for ax in (ax1, ax2):
    ax.axvline(TITIN, color=BLACK, lw=0.7, ls=":")
    ax.axvline(PKZ, color=BLACK, lw=0.7, ls="--")
ax1.text(TITIN + 400, 1500, "titin\n35,991 aa", fontsize=6.5, va="top")
ax1.text(PKZ + 400, 1500, "PKZILLA-1\n45,212 aa", fontsize=6.5, va="top")
lab = {"KAL1508129.1": "PKZILLA-1", "UPI002B62E083": "Omnitrophota\nORF", "KAI4551101.1": "KAI4551101.1", "KAK2904762.1": "KAK2904762.1",
       "WP_406270362.1": "Nocardia\nNBC_00881", "WP_327112628.1": "Nocardia NBC_01730", "XP_065194412.1": "Sycon"}
for r in main:
    if r["rep_accession"] in lab:
        x, y = int(r["length"]), rows[r["grade"]]
        ax2.annotate(lab[r["rep_accession"]], (x, y), xytext=(x, y + 0.45), fontsize=6, ha="center",
                     arrowprops=dict(arrowstyle="-", lw=0.4, color=GREY))
ax2.set_xlim(35000, 100000)
ax2.set_xlabel("Protein length (aa)")
ax2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{int(v):,}"))
panel(ax2, "b")
save(fig, "fig1_lengths_grades")

# ---- Figure 2: grade composition by evidence stage -----------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.2))
stages = [("grade_stage1", "1  sequence + locus"), ("grade_stage2_rna", "2  + RNA-seq junctions"),
          ("grade", "3  + source-read check")]
for ax, (title, rows_) in zip(axes, (("All 510 sequences > 35,991 aa", main),
                                     ("57 sequences > 45,212 aa", [r for r in main if int(r["length"]) > PKZ]))):
    n = len(rows_)
    for i, (col, lab_) in enumerate(stages):
        c = Counter(r[col] for r in rows_)
        left = 0
        for gr in "ABCD":
            w = c[gr] / n
            ax.barh(2 - i, w, left=left, color=GC[gr], edgecolor="white", linewidth=0.6, height=0.7)
            if c[gr]:
                txt = str(c[gr])
                if w > 0.035:
                    ax.text(left + w / 2, 2 - i, txt, ha="center", va="center", fontsize=6.5,
                            color="white" if gr in "AD" else BLACK)
                else:
                    ax.text(left + w / 2, 2 - i + 0.42, txt, ha="center", va="bottom", fontsize=6.5, color=BLACK)
            left += w
    ax.set_yticks([2, 1, 0])
    ax.set_yticklabels([s[1] for s in stages])
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(["0", "25", "50", "75", "100 %"])
    ax.set_title(title, loc="left")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
axes[1].set_yticklabels([])
handles = [matplotlib.patches.Patch(color=GC[g], label=GRADE_NAME[g]) for g in "ABCD"]
fig.legend(handles=handles, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.12))
panel(axes[0], "a")
panel(axes[1], "b")
fig.tight_layout()
save(fig, "fig2_grades_by_stage")

# ---- Figure 3: a translated tandem repeat --------------------------------------------------
prot = fasta_gz(os.path.join(D, "fasta", "unique_ge35992.fasta.gz"))
md5 = {r["rep_accession"]: r["md5"] for r in lc}


def self_identity(s, lags):
    a = np.frombuffer(s.encode(), dtype=np.uint8)
    return [float((a[:-k] == a[k:]).mean()) for k in lags]


fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.3))
lags = list(range(1, 121))
for acc, lab_, col in (("KAI4551101.1", "KAI4551101.1 (98,182 aa)", GC["D"]),
                       ("UPI002B62E083", "UPI002B62E083 (85,803 aa)", BLUE),
                       ("KAL1508129.1", "PKZILLA-1 (45,212 aa)", GREY)):
    axes[0].plot(lags, self_identity(prot[md5[acc]], lags), color=col, lw=0.9, label=lab_)
axes[0].set_xlabel("Lag (aa)")
axes[0].set_ylabel("Identity of the protein with itself")
axes[0].set_ylim(0, 1)
axes[0].set_title("Protein self-identity by lag", loc="left")
panel(axes[0], "a")
# DNA: the annotated CDS and its stops in all six frames
locus = fasta_gz(os.path.join(D, "manual", "KAI4551101_locus.fa.gz"))
lseq = next(iter(locus.values()))
# full CDS join from the GenPept record (the parsed table keeps a truncated copy)
gp = gzip.open(os.path.join(D, "raw", "ncbi_gp.txt.gz"), "rt").read()
i0 = gp.find("/coded_by=", gp.find("KAI4551101"))
coded_by = re.sub(r"\s+", "", gp[i0 + 11:gp.find('"', i0 + 11)])
parts = [(int(a), int(b)) for a, b in re.findall(r":<?(\d+)\.\.>?(\d+)", coded_by)]
off = int(next(iter(locus)).split(":")[1].split("-")[0])
cds = "".join(lseq[a - off:b - off + 1] for a, b in sorted(parts))
if coded_by.startswith("complement"):
    cds = cds.translate(COMP)[::-1]
assert len(cds) == 3 * 98182 + 3, len(cds)
mags = fasta_gz(os.path.join(D, "manual", "DAOTDT01.fasta.gz"), key=lambda h: h.split()[0].split("|")[-1])
tdt = json.load(open(os.path.join(D, "dna_checks", "omnitrophota_DAOTDT_HRZ03494.json")))
orf = mags[tdt["contig"]][tdt["gene"][0] - 1:tdt["gene"][1]]


def stops_per_k(s):
    out = []
    for strand in (s, s.translate(COMP)[::-1]):
        for f in range(3):
            cod = [strand[i:i + 3] for i in range(f, len(strand) - 2, 3)]
            out.append(1000 * sum(c in ("TAA", "TAG", "TGA") for c in cod) / len(cod))
    return out


x = np.arange(6)
w = 0.38
for i, (seq, lab_, col) in enumerate(((cds, "KAI4551101.1 CDS", GC["D"]), (orf, "UPI002B62E083 ORF", BLUE))):
    v = stops_per_k(seq)
    axes[1].bar(x + (i - 0.5) * w, v, width=w, color=col, label=lab_)
    for xi, vi in zip(x, v):
        if vi < 3:
            axes[1].text(xi + (i - 0.5) * w, vi + 1, f"{vi:.1f}", ha="center", fontsize=5.5, rotation=90)
axes[1].set_xticks(x)
axes[1].set_xticklabels(["+1\n(annot.)", "+2", "+3", "−1", "−2", "−3"])
axes[1].set_ylabel("Stop codons per 1,000 codons")
axes[1].set_title("Stops in the six frames", loc="left")
panel(axes[1], "b")
dl = list(range(1, 61))
for seq, lab_, col in ((cds, "KAI4551101.1 CDS", GC["D"]), (orf, "UPI002B62E083 ORF", BLUE)):
    axes[2].plot(dl, self_identity(seq[:200000], dl), color=col, lw=0.9, label=lab_)
axes[2].axhline(0.25, color=GREY, lw=0.5, ls=":")
axes[2].set_xlabel("Lag (bp)")
axes[2].set_ylabel("Identity of the DNA with itself")
axes[2].set_ylim(0, 1)
axes[2].set_title("CDS self-identity by lag", loc="left")
panel(axes[2], "c")
h3 = [matplotlib.lines.Line2D([], [], color=GC["D"], lw=2, label="KAI4551101.1 (98,182 aa; grade D)"),
      matplotlib.lines.Line2D([], [], color=BLUE, lw=2, label="UPI002B62E083 (85,803 aa; grade B)"),
      matplotlib.lines.Line2D([], [], color=GREY, lw=2, label="PKZILLA-1 (45,212 aa; panel a only)")]
fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.legend(handles=h3, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.0))
save(fig, "fig3_tandem_repeat")

# ---- Figure 4: the merged fish titin pair -------------------------------------------------
fp = [r for r in tsv("followups/fish_titin_pair.tsv") if r["model"] == "KAK2904762.1"]
fig = plt.figure(figsize=(7.0, 4.6))
gs4 = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1.7, 1])
a1 = fig.add_subplot(gs4[0, :])
a2 = fig.add_subplot(gs4[1, 0], sharex=a1)
a3 = fig.add_subplot(gs4[1, 1])
fmt = matplotlib.ticker.FuncFormatter(lambda v, p: f"{int(v):,}")
for par, col, mk in (("ttna", BLUE, "o"), ("ttnb", ORANGE, "s")):
    xs = [int(r["window_start_aa"]) + 500 for r in fp if r["best_zebrafish_paralogue"] == par]
    ys = [int(r[par + "_pos"].split("-")[0]) for r in fp if r["best_zebrafish_paralogue"] == par]
    a1.scatter(xs, ys, s=12, color=col, marker=mk, edgecolor="none", label=f"best match: zebrafish {par}")
a1.set_ylabel("Matched position in\nzebrafish titin (aa)")
a1.yaxis.set_major_formatter(fmt)
a1.legend(loc="upper left")
a1.set_title("KAK2904762.1 (Cirrhinus molitorella, 59,714 aa): 1,000-aa windows searched against zebrafish ttna and ttnb",
             loc="left")
panel(a1, "a")
tot, exons = {}, None
with tarfile.open(os.path.join(D, "rnaseq_run_results.tar.gz")) as t:
    for m in t.getmembers():
        if m.name.startswith("results/Cirrhinus_molitorella") and m.name.endswith(".json") \
                and not os.path.basename(m.name).startswith("._"):
            d = json.load(t.extractfile(m))
            for L in d["loci"]:
                if L["rep_accession"] == "KAK2904762.1":
                    exons = exons or L["exons"]
                    for i in L["introns"]:
                        tot[(i[0], i[1])] = tot.get((i[0], i[1]), 0) + i[2]
strand = next(r["strand"] for r in lc if r["rep_accession"] == "KAK2904762.1")
ex = sorted((e[0], e[1]) for e in exons)
ex = ex[::-1] if strand == "-" else ex
cum, pts = 0, []
for i, (a, b) in enumerate(ex[:-1]):
    cum += b - a + 1
    nxt = ex[i + 1]
    key = (b + 1, nxt[0] - 1) if strand == "+" else (nxt[1] + 1, a - 1)
    pts.append((cum / 3, tot.get(key, 0)))
px = np.array([p[0] for p in pts])
py = np.array([p[1] for p in pts])
j = int(np.argmin(np.abs(px - 30675)))
a2.scatter(px, py + 1, s=4, color=BLACK, edgecolor="none")
a2.scatter([px[j]], [py[j] + 1], s=30, facecolor="none", edgecolor=GC["D"], lw=1.0, zorder=5)
a2.set_yscale("log")
a2.yaxis.set_major_formatter(fmt)
a2.set_ylabel("Junction reads + 1\n(each intron; unique, 6 runs)")
a2.set_xlabel("Position in the model (aa)")
a2.xaxis.set_major_formatter(fmt)
panel(a2, "b")
for ax in (a1, a2):
    ax.axvline(31001, color=BLACK, lw=0.7, ls="--")
a1.text(31400, a1.get_ylim()[1] * 0.97, "second titin copy starts (aa 31,001)", fontsize=6.5, va="top")
m = (px > 29000) & (px < 33000)
a3.vlines(px[m], 0.8, py[m] + 1, color=GREY, lw=0.8)
a3.scatter(px[m], py[m] + 1, s=8, color=BLACK, edgecolor="none", zorder=4)
a3.scatter([px[j]], [py[j] + 1], s=40, facecolor="none", edgecolor=GC["D"], lw=1.1, zorder=5)
a3.set_yscale("log")
a3.set_ylim(0.8, 600)
a3.yaxis.set_major_formatter(fmt)
a3.axvline(31001, color=BLACK, lw=0.7, ls="--")
a3.annotate(f"intron at aa {int(px[j]):,}:\n{py[j]} junction reads", (px[j], 1.0), xytext=(29150, 60), fontsize=6.5,
            arrowprops=dict(arrowstyle="-", lw=0.5, color=GC["D"]))
a3.set_xlim(29000, 33000)
a3.xaxis.set_major_formatter(fmt)
a3.set_xlabel("Position in the model (aa)")
a3.set_title("Zoom on the join", loc="left")
panel(a3, "c")
fig.tight_layout()
save(fig, "fig4_fish_titin_pair")

# ---- Figure 5: the 85,803-aa Omnitrophota ORF ---------------------------------------------
dna = {k: json.load(open(os.path.join(D, "dna_checks", f"omnitrophota_{k}.json")))
       for k in ("DAOTDT_HRZ03494", "DAOPMG_HPD84636", "DAOOUB_HOW35280")}
labels = {"DAOTDT_HRZ03494": "DAOTDT01 (SRR5007146)", "DAOPMG_HPD84636": "DAOPMG01 (SRR5007150)",
          "DAOOUB_HOW35280": "DAOOUB01 (SRR4244739)"}
cols = {"DAOTDT_HRZ03494": BLUE, "DAOPMG_HPD84636": ORANGE, "DAOOUB_HOW35280": PURPLE}


def to_cds(j):
    gs, ge = j["gene"]
    n = ge - gs + 1
    il = j["illumina"]
    nw = len(il["window_depth_unique"])
    wc = np.array([1000 * i + 500 for i in range(nw)], float)
    bd = np.array([1000 * (i + 1) for i in range(nw - 1)], float)
    if j["strand"] == "-":
        wc, bd = n - wc, n - bd
    return wc / 1000, bd / 1000, il


num = json.load(open(os.path.join(ROOT, "paper", "numbers.json")))
ra, rb = num["rep_copy_a"], num["rep_copy_b"]
fig = plt.figure(figsize=(7.0, 5.6))
gsp = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.25], width_ratios=[1.6, 1])
ax_d = fig.add_subplot(gsp[0, :])
ax_s = fig.add_subplot(gsp[1, :], sharex=ax_d)
ax_z = fig.add_subplot(gsp[2, 0])
ax_p = fig.add_subplot(gsp[2, 1])
wc0, _, il0 = to_cds(dna["DAOTDT_HRZ03494"])
ax_d.plot(wc0, il0["window_depth_all_mapq"], color=BLACK, lw=0.7, label="all alignments")
ax_d.plot(wc0, il0["window_depth_unique"], color=BLUE, lw=1.0, label="unique alignments (MAPQ ≥ 20)")
for k, j in dna.items():
    wc, bd, il = to_cds(j)
    o = np.argsort(bd)
    ax_s.plot(bd[o], np.array(il["spanning_per_boundary"])[o], color=cols[k], lw=0.8, label=labels[k])
    E = np.array([e if e is not None else np.nan for e in il["expected_spanning_per_boundary"]])
    ax_s.plot(bd[o], E[o], color=cols[k], lw=0.6, ls=":")
for ax in (ax_d, ax_s):
    for (a, b) in (ra, rb):
        ax.axvspan(a / 1000, b / 1000, color=GREY, alpha=0.35, lw=0)
ax_d.set_ylabel("Read depth\n(1-kb windows)")
ax_d.set_ylim(0, 75)
ax_d.legend(loc="lower left", ncol=2, fontsize=6, bbox_to_anchor=(0.0, 0.02))
ax_d.set_title("Read depth, DAOTDT01 (SRR5007146); grey bands: the two copies of the internal repeat",
               loc="left", fontsize=7)
panel(ax_d, "a")
ax_s.set_ylabel("Unique pairs spanning\neach window boundary")
ax_s.set_ylim(-12, 80)
ax_s.set_xlabel("Position in the coding sequence (kb from the start codon)")
ax_s.set_title("Solid: observed in each MAG's own reads; dotted: expected from simulated library pairs", loc="left",
               fontsize=7)
ax_s.legend(loc="lower left", ncol=3, fontsize=6, bbox_to_anchor=(0.0, 0.0))
ax_s.set_xlim(0, 258)
panel(ax_s, "b")
for k, j in dna.items():
    wc, bd, il = to_cds(j)
    o = np.argsort(bd)
    m = (bd[o] > 215) & (bd[o] < 243)
    ax_z.plot(bd[o][m], np.array(il["spanning_per_boundary"])[o][m], color=cols[k], lw=0.9, marker="o", ms=2.5)
    E = np.array([e if e is not None else np.nan for e in il["expected_spanning_per_boundary"]])
    ax_z.plot(bd[o][m], E[o][m], color=cols[k], lw=0.7, ls=":")
    for b in il["boundaries_untestable_for_library"]:
        x0 = (b - j["gene"][0]) if j["strand"] == "+" else (j["gene"][1] - b + 1)
        ax_z.scatter([x0 / 1000], [0], s=40, marker="x", color=cols[k], zorder=5, lw=1.1)
for (a, b), t in ((ra, "copy A"), (rb, "copy B")):
    ax_z.axvspan(a / 1000, b / 1000, color=GREY, alpha=0.35, lw=0)
    ax_z.text((a + b) / 2000, ax_z.get_ylim()[1] if False else 62, t, ha="center", fontsize=6.5)
ax_z.set_xlim(215, 243)
ax_z.set_ylim(-4, 68)
ax_z.set_xlabel("CDS position (kb)")
ax_z.set_ylabel("Unique spanning pairs")
ax_z.set_title("Around the repeat; ×: boundary inside a copy (expected 0, observed 0)", loc="left", fontsize=7)
panel(ax_z, "c")
# dot plot of the repeat region against itself (exact 30-mers)
lo, hi = 219000, 239000
reg = orf[lo:hi]
K = 30
idx = {}
for i in range(len(reg) - K + 1):
    idx.setdefault(reg[i:i + K], []).append(i)
xs, ys = [], []
for v in idx.values():
    if len(v) > 1:
        for a in v:
            for b in v:
                if a != b:
                    xs.append((lo + a) / 1000)
                    ys.append((lo + b) / 1000)
ax_p.plot([lo / 1000, hi / 1000], [lo / 1000, hi / 1000], color=GREY, lw=0.6)
ax_p.scatter(xs, ys, s=0.6, color=BLACK, edgecolor="none")
ax_p.set_xlim(lo / 1000, hi / 1000)
ax_p.set_ylim(lo / 1000, hi / 1000)
ax_p.set_aspect("equal")
ax_p.set_xlabel("CDS position (kb)")
ax_p.set_ylabel("CDS position (kb)")
ax_p.set_title(f"Shared 30-mers (diagonal = self); repeat {num['rep_lag_bp']:,} bp apart", loc="left", fontsize=7)
panel(ax_p, "d")
fig.tight_layout()
save(fig, "fig5_omnitrophota_orf")

# ---- Figure 6: controls of the grading rules -------------------------------------------------
CT = os.path.join(ROOT, "controls")
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [1, 1.15]})
ctl = list(csv.DictReader(open(os.path.join(CT, "d1", "d1_controls.tsv")), delimiter="\t"))
hold = list(csv.DictReader(open(os.path.join(CT, "d1", "holdout.tsv")), delimiter="\t"))
ax_a.scatter([float(r["rep25"]) for r in ctl], [float(r["period_identity"]) for r in ctl], s=7, color=GREY, alpha=0.6,
             edgecolor="none", label=f"real proteins, design set (n = {len(ctl)})")
real_fl = [r for r in ctl if r["D1"] == "True"] + [r for r in hold if r["D1"] == "True"]
ax_a.scatter([float(r["rep25"]) for r in real_fl], [float(r["period_identity"]) for r in real_fl], s=16, color=BLUE,
             marker="s", edgecolor="none", label=f"real proteins with the pattern (n = {len(real_fl)})")
dd = [r for r in main if "D1 repeat translation" in r["grade_reason_stage1"]]
ff = [r for r in main if "repeat flag" in r["grade_reason_stage1"]]
ax_a.scatter([float(r["rep25_frac"]) for r in dd], [float(r["protein_period_identity"]) for r in dd], s=14, color=GC["D"],
             edgecolor="none", label=f"audit, pattern + DNA contradiction: D (n = {len(dd)})")
ax_a.scatter([float(r["rep25_frac"]) for r in ff], [float(r["protein_period_identity"]) for r in ff], s=14,
             facecolor="none", edgecolor=GC["C"], lw=0.8, label=f"audit, pattern only: C, flagged (n = {len(ff)})")
ax_a.plot([0, 0.8, 0.8], [0.6, 0.6, 0], color=BLACK, lw=0.6, ls="--")
ax_a.text(0.02, 0.56, "below and left of the dashed line:\nno repeat pattern", fontsize=6, va="top")
ax_a.set_xlim(0, 1.02)
ax_a.set_ylim(0, 1.02)
ax_a.set_xlabel("Share of residues in repeated 25-mers (rep25)")
ax_a.set_ylabel("Self-identity at the best period")
ax_a.legend(loc="upper left", fontsize=5.6, handletextpad=0.3, bbox_to_anchor=(0.0, 1.0))
ax_a.set_title("Repeat pattern in real proteins and in the audit", loc="left")
panel(ax_a, "a")
pw = list(csv.DictReader(open(os.path.join(CT, "rna_power", "rna_power.tsv")), delimiter="\t"))
fr = sorted({float(r["fraction_wrong"]) for r in pw})
for kind, col in (("shift", BLUE), ("skip", ORANGE)):
    frk = [f for f in fr if not (kind == "skip" and f == 1.0)]     # all exons skipped leaves one intron: untestable
    for o, ls in (("v1", ":"), ("v2", "-")):
        y = []
        for f in frk:
            rr = [r for r in pw if r["error"] == kind and float(r["fraction_wrong"]) == f]
            tot = sum(int(r["replicates"]) for r in rr)
            y.append(sum(int(r[f"{o}_contradicted"]) + int(r[f"{o}_weak, under-supported where covered"]) for r in rr) / tot)
        ax_b.plot(frk, y, color=col, ls=ls, lw=1.1, marker="o" if o == "v2" else None, ms=2.5,
                  label=f"{'introns shifted 15-60 bp' if kind == 'shift' else 'exons skipped'}, rules {o}")
    yc = []
    for f in frk:
        rr = [r for r in pw if r["error"] == kind and float(r["fraction_wrong"]) == f]
        tot = sum(int(r["replicates"]) for r in rr)
        yc.append(sum(int(r["v2_contradicted"]) for r in rr) / tot)
    ax_b.plot(frk, yc, color=col, lw=0.8, ls="--", label=f"{'shifted' if kind == 'shift' else 'skipped'}: 'contradicted' only")
ax_b.set_xlabel("Share of the model's introns made wrong")
ax_b.set_ylabel("Share of wrong models downgraded")
ax_b.set_ylim(-0.02, 1.02)
ax_b.legend(loc="lower right", fontsize=5.6)
ax_b.set_title(f"RNA-seq rule on wrong models ({len({r['model'] for r in pw})} real models)", loc="left")
panel(ax_b, "b")
fig.tight_layout()
save(fig, "fig6_controls")


# ---- Figure 7: peptides of the 39,686-aa Velamenicoccus protein (PXD025008 re-searches) -----------
pp = os.path.join(ROOT, "proteomics", "results", "raw60_giant_peptides.tsv")
if os.path.exists(pp):
    gseq = prot[md5["WP_128698900.1"]]

    def upeps(name):
        return [r for r in csv.DictReader(open(os.path.join(ROOT, "proteomics", "results", name)), delimiter="\t")
                if r["unique_to_giant"] == "True"]
    gpep, mpep = upeps("raw60_giant_peptides.tsv"), upeps("full6_giant_peptides.tsv")
    fig, (t1, t2) = plt.subplots(2, 1, figsize=(7.0, 2.9), sharex=True, gridspec_kw={"height_ratios": [1.4, 1]})
    for r in gpep:
        for a in r["start_positions"].split(";"):
            t1.vlines(int(a), 0, int(r["psms"]), color=BLUE, lw=0.8)
    raw_set = {r["peptide"] for r in gpep}
    for r in mpep:
        if r["peptide"] not in raw_set:
            a = int(r["start_positions"].split(";")[0])
            t1.plot([a], [int(r["psms"])], marker="o", ms=3, mfc="none", mec=ORANGE, mew=0.8, ls="none",
                    label="only in the peak-list search" if r is mpep[0] or True else None)
    h, l_ = t1.get_legend_handles_labels()
    if h:
        t1.legend(h[:1], l_[:1], loc="upper right", fontsize=6, frameon=False)
    t1.set_ylabel("PSMs per\nunique peptide")
    t1.set_title(f"WP_128698900.1 (39,686 aa): {len(gpep)} unique peptides at 1 % FDR in the 60 raw files of PRIDE PXD025008",
                 loc="left")
    vx = json.load(open(os.path.join(ROOT, "proteomics", "results", "verify_extra.json")))
    blocks_ = [(a - 1, b - 1) for a, b in vx["repeat_blocks_ge100aa"]]  # 20-mers occurring twice, blocks >= 100 aa
    for ax in (t1, t2):
        for a, b in blocks_:
            ax.axvspan(a + 1, b + 1, color=ORANGE, alpha=0.25, lw=0, zorder=0)
    edges = np.arange(0, len(gseq) + 1000, 1000)
    th = np.array(vx["theoretical_tryptic_7_40_per_1000aa"])
    ob, _ = np.histogram([int(a) for r in gpep for a in r["start_positions"].split(";")], bins=edges)
    t2.bar(edges[:-1], th, width=1000, align="edge", color=GREY, alpha=0.5, edgecolor="white", lw=0.3, zorder=1,
           label="fully tryptic peptides of 7-40 aa in the sequence")
    t2.bar(edges[:-1], ob, width=1000, align="edge", color=BLUE, edgecolor="white", lw=0.3, zorder=2,
           label="unique peptides identified (raw files)")
    t2.set_ylim(0, 90)
    t2.set_ylabel("Peptides per\n1,000 aa")
    t2.set_xlabel("Position in WP_128698900.1 (aa); shaded: internal repeat blocks")
    t2.legend(loc="upper right", fontsize=6, ncol=2)
    t2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{int(v):,}"))
    t2.set_xlim(0, len(gseq))
    t2.set_xticks(list(range(0, len(gseq), 5000)))
    panel(t1, "a")
    panel(t2, "b")
    fig.tight_layout()
    save(fig, "fig7_peptides")

print("figures written to", OUT)
