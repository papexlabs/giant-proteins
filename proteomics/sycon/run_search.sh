#!/bin/bash
# PXD060105 (Sycon ciliatum spicule matrix) re-search: RAW -> indexed mzML (ThermoRawFileParser, vendor
# peak picking), database (build_db.py), Sage with the authors' stated settings where given (trypsin,
# 1 missed cleavage, 10 ppm precursor, carbamidomethyl C fixed, Met oxidation variable) and a 20 ppm
# fragment tolerance for the Orbitrap HCD spectra; then analyze.py. 4 threads, telemetry off.
set -eu
cd "$(dirname "$0")"
mkdir -p mzml search/out results
for f in raw/*.raw; do b=$(basename "$f" .raw); [ -s "mzml/$b.mzML" ] || ../tools/trfp/ThermoRawFileParser -i "$f" -o mzml -f 2 -m 0 -l 1 >> mzml/convert.log 2>&1; done
python3 -I build_db.py
export RAYON_NUM_THREADS=4
/usr/bin/time -l ../tools/sage-v0.14.7-aarch64-apple-darwin/sage --disable-telemetry-i-dont-want-to-improve-sage \
  --batch-size 1 search/sage_config.json mzml/*.mzML 2>&1 | tee search/sage.log
python3 -I analyze.py
