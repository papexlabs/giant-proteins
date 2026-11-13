#!/bin/bash
# Download the 60 Thermo RAW files of PXD025008 into raw_thermo/ (resumable, 4 in parallel).
# File names, sizes and SHA-1 checksums come from results/pride_files_PXD025008.json (PRIDE API);
# a file counts as complete only when its size and SHA-1 match. Re-run to resume.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/raw_thermo"
export UA='PapexResearchBot/1.0 (+https://papex.org)'
export BASE=https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD025008
python3 -I -c "
import json
for x in json.load(open('$HERE/results/pride_files_PXD025008.json')):
    if (x['fileCategory'] or {}).get('value') == 'RAW':
        print(x['fileName'], x['fileSizeBytes'], x['checksum'])" > want.txt
one() {
  f=$1; want=$2; sha=$3
  for i in $(seq 1 40); do
    have=$(stat -f %z "$f" 2>/dev/null || stat -c %s "$f" 2>/dev/null || echo 0)
    if [ "$have" = "$want" ]; then
      [ "$(shasum -a 1 "$f" | cut -d' ' -f1)" = "$sha" ] && { echo "OK $f"; return 0; }
      echo "SHA1 MISMATCH $f, refetching"; rm -f "$f"; have=0
    fi
    [ "$have" -gt "$want" ] && rm -f "$f"
    curl -sSf -C - --connect-timeout 30 --speed-limit 50000 --speed-time 120 -A "$UA" -o "$f" "$BASE/$f" || sleep $((i * 5))
  done
  echo "FAILED $f"
}
export -f one
xargs -P 4 -L 1 bash -c 'one "$0" "$1" "$2"' < want.txt > download.log 2>&1
grep -c '^OK' download.log
