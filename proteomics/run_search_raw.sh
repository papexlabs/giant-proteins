#!/bin/bash
# Independent search of all 60 PXD025008 RAW files (Extra 2, 13 Nov 2026).
# 1. Waits until run_download_raw.sh has verified all 60 files (size + SHA-1).
# 2. Converts each RAW file to indexed mzML with ThermoRawFileParser (tools/trfp, release
#    v.2.0.0-dev for macOS arm64, unpacked in the project and ad-hoc code-signed locally so macOS
#    runs it; vendor peak picking on, the converter's default), metadata JSON alongside.
# 3. Searches them with the SAME Sage binary, database (db/search_db.fasta) and settings as the
#    peak-list search (search/sage_config.json); only the input files and output folder differ.
set -eu
cd "$(dirname "$0")"
n=$(wc -l < raw_thermo/want.txt | tr -d ' ')
until [ "$(grep -c '^OK' raw_thermo/download.log 2>/dev/null || echo 0)" -ge "$n" ]; do
  grep -q '^FAILED' raw_thermo/download.log && { echo "a RAW download FAILED; stopping"; exit 1; }
  sleep 60
done
mkdir -p mzml search/out_raw
for f in raw_thermo/*.raw; do
  b=$(basename "$f" .raw)
  [ -s "mzml/$b.mzML" ] || tools/trfp/ThermoRawFileParser -i "$f" -o mzml -f 2 -m 0 -l 1 >> mzml/convert.log 2>&1
done
python3 -I - <<'PY'
import json
c = json.load(open("search/sage_config.json"))
c["output_directory"] = "search/out_raw"
json.dump(c, open("search/sage_config_raw.json", "w"), indent=1)
PY
export RAYON_NUM_THREADS=4
/usr/bin/time -l tools/sage-v0.14.7-aarch64-apple-darwin/sage --disable-telemetry-i-dont-want-to-improve-sage \
  --batch-size 1 search/sage_config_raw.json mzml/*.mzML 2>&1 | tee search/sage_raw.log
python3 -I analyze_verified.py search/out_raw raw60
