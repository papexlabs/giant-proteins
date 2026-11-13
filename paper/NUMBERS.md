# Where each number in the manuscript comes from

The LaTeX sources (`paper/manuscript.tex`, `paper/supplement.tex`) print every number taken from the data
through `\nm{key}`. The values come from `paper/numbers.tex`, which `paper/make_numbers.py` writes from the
repository's files on every build; an unknown key stops the build. For any number in the PDF:

1. find its `\nm{key}` in the `.tex` source;
2. look the key up in `paper/NUMBERS_generated.tsv` (columns key, value, source), which names the file and
   column or the computation it comes from.

`paper/numbers.json` holds the same values in machine-readable form. Printing rules (thousands separators,
q-values to two significant digits, ranges with an en dash) are in `tex_value()` in `make_numbers.py`.

## Keys whose source is not a data file

| Key | Value | Source |
|---|---|---|
| `s1_fix_B_to_C` | 48 | [notes] `notes/FEASIBILITY.md` §14, `notes/RESULTS.md` §2 (Stage 1 counts before and after the isoform fix) |
| `titin_downgrades_withdrawn` | 5 | [notes] `notes/RESULTS.md` §3 |
| `stage1_cpu_hours` | 4.4 | [notes] `notes/FEASIBILITY.md` §11 |
| `raw_wall_min` | 56 | `proteomics/search/sage_raw.log` ("finished in 3337s"; the log is not committed, the value is) |
| `prot_peak_last_q` | 0.0096 | `proteomics/results/full6_giant_peptides.tsv`, peptide LNTIVTQNK |
| derived keys (`orf_len_aa`, `orf_len_if_excised`, `raw_giant_unc_pct`, `cirrhinus_first_end`, `rep_copy_bp_kb`, `pw_*`) | | arithmetic on measured keys, stated in the source column of `NUMBERS_generated.tsv` |

## Numbers written directly in the text

These are not measurements of this study, so they are not keys.

| Text | Source |
|---|---|
| 7 of 17 PKZILLA-1 exons with unique peptides | Fallon et al. 2024, full text (PMC11416037, checked 13 Nov 2026) |
| 85,804 aa (West-Roberts et al.) | bioRxiv abstract, doi:10.1101/2023.11.21.568195 (checked 13 Nov 2026) |
| more than 41,000 UniProtKB proteins of at least 5,000 aa | Amaral & Devos 2024, abstract |
| "4,384 kDa", "39,678 amino acids" | Kizina et al. 2022, full text (`proteomics/paper_text/PMC9004380_full.txt`) |
| 16,000 ("beyond residue 16,000"), 6,000 ("first 6,000 aa") | round positions chosen for the text; the measured positions are keys |
| rule thresholds (0.6, 0.8, 5 %, 50 %, 2 frames, 0.95, 0.74, 70 %, ±20 %, ±5 %, 95 %, 90 %, 3 reads, MAPQ 10 / 20, depth 10, 50 reads, 20 %, 1 kb, 5 pairs, 3 long reads, 500-bp flanks, 99 %, 10 %, 0.5–2, √2, 3 kb, 60 million reads, 1,800 s, 0.2–6 GB) | `scripts/13_grade.py`, `cloud/summarize.py`, `scripts/15_apply_rnaseq.py`, `cloud/dna_check.py` (`VERDICT`), `cloud/driver.py`, `cloud/make_jobs.py` |
| junction 50 bp either side; 10-mer share ≥ 50 % for Sycon fragments; fifths of the length | `controls/rna_rep/junction_mappability.py`, `proteomics/sycon/build_db.py`, `scripts/13_grade.py` (S6) |
| search settings (trypsin, 1 missed cleavage, 7–40 aa, 10 ppm, 0.5 Da, 20 ppm, isotope errors 0–1, ≤ 2 Met oxidations, 4 threads) | `proteomics/search/sage_config.json`, `proteomics/sycon/search/sage_config.json`, `proteomics/run_search*.sh` |
| control construction (10-kb deletion, 4.55-kb collapse, 5-kb duplication, 13.8-kb segmental repeat, 9.37 kb, 9.0-kb repeat, 1.1 kb, 15–60 bp shifts, 100 testable introns, 1,000–4,999 aa hold-out, 300 / 2,000 / 2,500 / 5,000 aa design-set limits) | `controls/dna/build_refs.py`, `controls/dna/summary.json`, `controls/rna_power/run_rna_power.py`, `controls/d1/run_d1_controls.py`, `run_d1_holdout.py` |
| 1.2-kb and 1.3 kb repeat (Omnitrophota), 10.4-kb segment | rounded from `rep_copy_bp`, the verdict text in `data/dna_checks/*.json`, and `rep_lag_bp` |
| 9 codons (QAT16263.1 starts 9 codons later) | `proteomics/db/WP_128698900.1.fasta` vs `QAT16263.1.fasta` (WP[10:] = QAT[1:]) |
| about 45 stops per 1,000 codons in the reverse frames | `ovis_stops_six_frames` (45.2–45.6) |
| software versions (minimap2 2.31, samtools 1.25, DIAMOND 2.2.8, PyHMMER 0.12.3, Sage 0.14.7, ThermoRawFileParser v.2.0.0-dev) | tool `--version` output on the Mac; `cloud/bootstrap.sh`; `proteomics/tools/` |
| dates (11 Nov 2026 snapshot, release 2026_04, 13 Nov file times) | `data/ncbi_counts.tsv`, `data/uniprotkb_ge35992.tsv` headers, `cloud/bundle_rep_rna/RULES.txt` |

## History

Until 13 Nov 2026 the manuscript was `paper/manuscript.md` with a hand-written number table here; it was
replaced by the LaTeX sources (see git history for the earlier version).
