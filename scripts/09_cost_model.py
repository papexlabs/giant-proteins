"""Cost model for Stages S2/S3 from the measured Stage 0 tables. All prices and
per-run CPU figures are ESTIMATES (parameters below); counts of species, loci,
SRA experiments and run sizes are MEASURED (data/*.tsv).

S2 (short-read RNA-seq junction support) applies to eukaryotic non-excluded loci
(T-arch / T-unclear) whose longest model lacks a RefSeq statement of all-intron
support. One pass per run baits reads for all loci of that species.
S3 (long reads) applies to candidate loci: all non-excluded loci > 45,212 aa plus
the longest titin-like locus per species (central), all non-excluded loci (high).
Output: printed table + data/cost_model.tsv
"""
import csv
import os
import statistics

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")
rows = list(csv.DictReader(open(os.path.join(D, "stage0_table.tsv")), delimiter="\t"))
avail = {r["organism"]: r for r in csv.DictReader(open(os.path.join(D, "species_data_availability.tsv")), delimiter="\t")}

loci = {}
for r in rows:
    if r["triage"] in ("T-arch", "T-unclear"):
        loci.setdefault(r["locus_key"], r)          # rows are sorted by length: first = longest
euk = [r for r in loci.values() if r["superkingdom"] == "Eukaryota"
       and r["rnaseq_samples_all_introns"] in ("", "0")]
s2_species = sorted({r["organism"] for r in euk})
with_rna = [s for s in s2_species if s in avail and int(avail[s]["sra_rnaseq_experiments"] or 0) > 0]
sizes = [float(avail[s]["median_rnaseq_run_size_gb"]) for s in with_rna if avail[s]["median_rnaseq_run_size_gb"]]
med_size = statistics.median(sizes) if sizes else 3.0

big = [r for r in loci.values() if int(r["length"]) > 45212]
cand_central = sorted({r["organism"] for r in big})
cand_high = sorted({r["organism"] for r in loci.values()})
lr = lambda sp: [s for s in sp if s in avail and int(avail[s]["sra_longread_genomic_experiments"] or 0) > 0]

P = {  # name: (low, central, high)  -- ESTIMATES
    "price_usd_per_vcpu_h_us_east_1_spot": (0.005, 0.008, 0.015),
    "s2_runs_per_species": (3, 5, 20),
    "s2_vcpu_h_per_run": (0.75, 1.5, 3.0),
    "s2_genome_index_vcpu_h_per_species": (0, 0, 16),
    "s3_vcpu_h_per_species": (10, 25, 60),
    "s3_gb_per_species": (20, 50, 150),
    "scratch_gb": (200, 500, 2000), "scratch_months": (0.5, 0.5, 1.0), "gp3_usd_gb_month": (0.09, 0.09, 0.09),
}
out = []
for i, case in enumerate(("low", "central", "high")):
    price = P["price_usd_per_vcpu_h_us_east_1_spot"][i]
    runs = sum(min(P["s2_runs_per_species"][i], int(avail[s]["sra_rnaseq_experiments"])) for s in with_rna)
    s2_vcpu = runs * P["s2_vcpu_h_per_run"][i] + len(with_rna) * P["s2_genome_index_vcpu_h_per_species"][i]
    s2_gb = runs * med_size
    s3_sp = lr(cand_central if case != "high" else cand_high)
    if case == "low":
        s3_sp = s3_sp[:10]
    s3_vcpu = len(s3_sp) * P["s3_vcpu_h_per_species"][i]
    s3_gb = len(s3_sp) * P["s3_gb_per_species"][i]
    other_vcpu = (20 + 50 + 50, 60 + 300 + 150, 200 + 1500 + 500)[i]   # S1b + S4 + S5 (see FEASIBILITY)
    storage = P["scratch_gb"][i] * P["scratch_months"][i] * P["gp3_usd_gb_month"][i]
    vcpu = s2_vcpu + s3_vcpu + other_vcpu
    out.append({"case": case, "s2_species_with_rnaseq": len(with_rna), "s2_runs": runs,
                "s2_gb": round(s2_gb), "s2_vcpu_h": round(s2_vcpu), "s3_species": len(s3_sp),
                "s3_gb": s3_gb, "s3_vcpu_h": s3_vcpu, "other_vcpu_h": other_vcpu,
                "total_vcpu_h": round(vcpu), "compute_usd": round(vcpu * price, 1),
                "storage_usd": round(storage, 1), "total_usd": round(vcpu * price + storage, 1)})
print(f"eukaryotic S2 loci {len(euk)} in {len(s2_species)} species; with >=1 SRA RNA-seq experiment: "
      f"{len(with_rna)}; median sampled run size {med_size:.2f} GB")
print(f"loci > 45,212 aa: {len(big)} in {len(cand_central)} organisms; with long-read genomic SRA data: "
      f"{len(lr(cand_central))}; all non-excluded organisms with long-read genomic data: {len(lr(cand_high))}/{len(cand_high)}")
cols = list(out[0].keys())
with open(os.path.join(D, "cost_model.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for o in out:
        f.write("\t".join(str(o[c]) for c in cols) + "\n")
        print(o)
