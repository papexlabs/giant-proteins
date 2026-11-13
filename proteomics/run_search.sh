#!/bin/bash
# Search the six PXD025008 MudPIT peak lists with Sage (open source) against db/search_db.fasta.
# Settings follow the authors' methods where stated (trypsin, 10 ppm precursor, ~1 Da fragment window for the
# LTQ ion trap, Met oxidation variable, no fixed modification); missed cleavages are limited to 1 and peptide
# length to 7-40 so the fragment index fits in this Mac's memory. 4 threads, telemetry off.
set -eu
cd "$(dirname "$0")"
python3 -I build_db.py
export RAYON_NUM_THREADS=4
/usr/bin/time -l tools/sage-v0.14.7-aarch64-apple-darwin/sage --disable-telemetry-i-dont-want-to-improve-sage \
  --batch-size 1 search/sage_config.json raw/Mudpit_*.MGF 2>&1 | tee search/sage.log
