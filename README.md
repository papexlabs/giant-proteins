# Evidence grades for protein database entries longer than titin

Code and data for an evidence-graded audit of every protein sequence longer than a 35,991-aa
human titin isoform in UniProtKB, UniParc and NCBI Protein (snapshots of 11 November 2026). Each sequence is
graded A–D by explicit checks, RNA-seq junction evidence, source-read checks of bacterial assemblies and
public proteomics. **Status: draft, not peer reviewed.**

The manuscript and its supplement (LaTeX sources and PDFs) will be added to `paper/` when the preprint is posted. The graded table
(Supplementary Table S1) is `paper/supplementary/Table_S1_leaderboard.tsv`, with its column dictionary beside it.

## Layout

| Folder | Contents |
|---|---|
| `scripts/` | Stages 0–1, grading, applying RNA-seq and DNA results, leaderboard (numbered in run order) |
| `cloud/` | RNA-seq and source-read checks that run on a cloud box: drivers, job bundles (`bundle*/`), launch and sync scripts |
| `controls/` | Controls of the grading rules: repeat rule (`d1/`), RNA-seq power test (`rna_power/`), known misassemblies (`dna/`), repeat-flag RNA stage (`rna_rep/`), threshold sensitivity (`sensitivity/`) |
| `proteomics/` | Re-search of PRIDE PXD025008 (peak lists and raw files) and PXD060105 (`sycon/`) |
| `data/` | Inventories, per-sequence checks and grades (`locus_checks.tsv`), RNA-seq results (`rnaseq_support.tsv`, all run results in `rnaseq_run_results.tar.gz`), DNA-check results (`dna_checks/`), leaderboard |
| `paper/` | LaTeX sources, figures, tables, numbers and the build |
| `release/` | Data manifest and release-tree scripts |

`DATA_MANIFEST.tsv` lists every stored data file with its size, every input the scripts download instead
(RNA-seq runs, genome assemblies, source reads, proteomics files, with sizes where the archive gives them) and
the git-ignored folders the scripts fill.

## Requirements

- Python 3.11+ (standard library for the pipeline); a virtual environment with `matplotlib`, `numpy` and
  `pyhmmer` (0.12.3) for figures and the domain scan: `python3 -m venv .venv && .venv/bin/pip install matplotlib numpy pyhmmer==0.12.3`
- minimap2 2.31, samtools (1.25 used), DIAMOND 2.2.8 for Stages 1–3
- Sage 0.14.7 and ThermoRawFileParser (v.2.0.0-dev used) for the proteomics re-searches; unpack them into
  `proteomics/tools/` (not redistributed here)
- tectonic for the PDFs
- Network access to UniProt, NCBI E-utilities and FTP, ENA and PRIDE. All requests use a generic
  User-Agent and stay within the providers' anonymous rate limits; no contact data is sent.

## Rerunning

Every stage reads the previous stage's files, so each can be rerun alone.

| Stage | Command(s) | Output |
|---|---|---|
| Inventory | `scripts/01_fetch_uniprot.py`, `02_fetch_ncbi.py`, `03_fetch_ncbi_records.py`, `04_flags.py`, `05_data_availability.py` | `data/*_ge35992.*`, `giant_unique_sequences.tsv`, `species_data_availability.tsv` |
| Stage 0 | `06_domain_scan.py`, `07_period_flags.py`, `08_stage0_table.py` | `stage0_table.tsv`, `stage0_loci.tsv` |
| Stage 1 | `10_locus_fetch.py`, `10b_controls.py`, `11_locus_dna_checks.py`, `12_orthology.py`, `13_grade.py` | `locus_checks.tsv` (grade_stage1) |
| Stage 2 (RNA-seq) | `cloud/make_jobs.py` builds a job bundle; `cloud/driver.py <bundle> <workdir>` runs it (on a cloud box via `cloud/launch.sh`, or locally with `--no-shutdown`); `cloud/summarize.py`; `scripts/15_apply_rnaseq.py` | `rnaseq_support.tsv`, `grade_stage2_rna` |
| Stage 3 (source reads) | `cloud/dna_check.py <job> <outdir>` per assembly (`cloud/bundle_dna/`); `scripts/20_apply_dna.py` | `data/dna_checks/*.json`, final `grade` |
| Leaderboard | `scripts/16_leaderboard.py` | `data/leaderboard.tsv` |
| Controls | `controls/d1/run_d1_controls.py`, `run_d1_holdout.py`; `controls/rna_power/run_rna_power.py`; `controls/dna/build_refs.py` then `run_controls.sh`; `controls/rna_rep/junction_mappability.py`; `controls/sensitivity/run_sensitivity.py` | `controls/*/` |
| Proteomics | `proteomics/run_download.sh`, `run_download_raw.sh`, `build_db.py`, `run_search.sh`, `run_search_raw.sh`, `analyze_verified.py`, `verify_extra.py`, `verify_raw.py`; `proteomics/sycon/run_download.sh`, `run_search.sh` | `proteomics/results/`, `proteomics/sycon/results/`, `data/followups/peptide_evidence.tsv` |
| Paper | `.venv/bin/python paper/build.py` (needs the LaTeX sources, added with the preprint) | `paper/manuscript.pdf`, `paper/supplement.pdf` |

To regrade from the stored results without any download, restore the synced RNA-seq results and run the
grading chain:

```bash
mkdir -p cloud_out && tar xzf data/rnaseq_run_results.tar.gz -C cloud_out
python3 scripts/13_grade.py && python3 cloud/summarize.py && python3 scripts/15_apply_rnaseq.py \
  && python3 scripts/20_apply_dna.py && python3 scripts/16_leaderboard.py
```

The threshold-sensitivity script checks that this chain reproduces the released grades exactly.

### Cloud runs

`cloud/launch.sh` starts one spot instance, pushes the code and a job bundle, and starts the driver; the box
arms its own hard stop at boot and terminates when it finishes. Account-specific settings (AWS CLI profile, key
pair, security group, subnets, AMI) are read from `cloud/aws.env`, which is git-ignored; copy
`cloud/aws.env.example` and fill in your own. `STAGE=rep cloud/launch.sh <HARD_STOP_UTC>` runs only the
repeat-flag RNA stage, under the rules fixed in `cloud/bundle_rep_rna/RULES.txt`.

## Numbers in the paper

Every data number in the PDFs is printed through `\nm{key}` from `paper/numbers.tex`, which
`paper/make_numbers.py` computes from the files above; `paper/NUMBERS_generated.tsv` gives each key's source and
`paper/NUMBERS.md` explains the mapping.

## Licence

Code (`scripts/`, `cloud/`, `release/` and every `.py` and `.sh` file) is under the MIT licence (`LICENSE`). Data, tables, figures and text are under CC BY 4.0 (`LICENSE-data.md`). Sequences and record fields copied from UniProt, NCBI, ENA and PRIDE keep their sources' terms.
