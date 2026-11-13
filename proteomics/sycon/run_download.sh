#!/bin/bash
# Download the 5 Thermo RAW files of PRIDE PXD060105 (Sycon ciliatum spicule matrix) into raw/.
# Names and sizes from raw/files.txt (PRIDE API file list, 13 Nov 2026; no checksums given):
# a file counts as complete when its size matches. Resumable.
set -u
cd "$(dirname "$0")/raw"
UA='PapexResearchBot/1.0 (+https://papex.org)'
grep '\.raw ' files.txt | while read -r f size loc; do
  url=$(echo "$loc" | sed 's#^ftp://#https://#')
  for i in $(seq 1 40); do
    have=$(stat -f %z "$f" 2>/dev/null || stat -c %s "$f" 2>/dev/null || echo 0)
    [ "$have" = "$size" ] && break
    [ "$have" -gt "$size" ] && rm -f "$f"
    curl -sSf -C - --connect-timeout 30 --speed-limit 50000 --speed-time 120 -A "$UA" -o "$f" "$url" || sleep $((i * 5))
  done
  echo "$f $(stat -f %z "$f" 2>/dev/null || stat -c %s "$f") of $size"
done
