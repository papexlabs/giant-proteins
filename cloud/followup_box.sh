#!/bin/bash
# Box-side runner for the follow-up stage (STAGE=followup in launch.sh):
#   1. DNA re-mapping checks (cloud/bundle_dna) -> $W/results/_dna/<job>.json
#   2. RNA-seq (cloud/bundle_followup_rna, existing RNA driver, saturation rule off) for the
#      three loci lost to the job-key collision (<species>_extra) and the four titin models
#      whose runs stopped on the saturation cap (<species>_nosat)
# The RNA driver writes "powering off" at the end and shuts the box down (terminate); the
# hard stop armed at boot stays in force either way. Both steps skip finished jobs/runs, so a
# relaunched box resumes. GP/W/RNA_BUNDLE/DNA_ARGS/RNA_ARGS can be overridden for a Mac test.
set -u
GP=${GP:-/home/ubuntu/gp}
W=${W:-/home/ubuntu/gpwork}
RNA_BUNDLE=${RNA_BUNDLE:-$GP/bundle_followup_rna}
mkdir -p "$W"
python3 "$GP/dna_driver.py" "$GP/bundle_dna" "$W" --par 3 --threads 32 ${DNA_ARGS:-} >> "$W/dna_driver.out" 2>&1
# --sat-cov 0: no saturation stop (the collision loci and the saturation-limited titin runs)
python3 "$GP/driver.py" "$RNA_BUNDLE" "$W" --threads 32 --par 8 --sat-cov 0 ${RNA_ARGS:-} >> "$W/driver.out" 2>&1
