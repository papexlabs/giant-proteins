#!/usr/bin/env bash
# Run cloud/dna_check.py on the control regions (jobs/ctl_*), Illumina only and Illumina + ONT.
# One mapping per reference: the first job on a reference maps the reads into work/<ref>/, later jobs
# reuse those BAMs (dna_check.py skips a mapping whose .ok marker exists). The in-silico chimeric
# control is skipped here (--no-negative-control): these regions are themselves the controls.
# Usage: controls/dna/run_controls.sh <out_subdir> [extra dna_check.py args]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$HERE/../.."
OUT="$HERE/results/${1:-v3}"
shift || true
mkdir -p "$OUT" "$HERE/work"
cd "$ROOT"                        # job files hold paths relative to the repository root
for ref in modified original; do
  for mode in ill illont; do
    for j in "$HERE"/jobs/ctl_*_"$mode".json; do
      if grep -q "\"assembly_path\": \".*/$ref.fa.gz\"" "$j"; then
        echo "$(date +%H:%M:%S) $(basename "$j" .json) on $ref"
        python3 "$ROOT/cloud/dna_check.py" "$j" "$OUT" --work "$HERE/work/$ref" --threads 4 --no-negative-control "$@" \
          2>&1 | grep -E "done|mapped|incomplete|Error|Traceback" || true
      fi
    done
  done
done
