#!/bin/bash
# Download the PXD025008 files used in the reanalysis (peak lists, authors' FASTA, SDRF, overview) into raw/.
# Resumable; a file counts as complete when its size equals the server's Content-Length.
set -u
cd "$(dirname "$0")/raw"
UA='PapexResearchBot/1.0 (+https://papex.org)'
BASE=https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD025008
FILES="
Mudpit_160211_OV3_P5_SK_JHa_373-Pe2_01.mzid_Mudpit_160211_OV3_P5_SK_JHa_373-Pe2_01.MGF
Mudpit_160212_OV3_P5_SK_JHa_373-Pe1_01.mzid_Mudpit_160212_OV3_P5_SK_JHa_373-Pe1_01.MGF
Mudpit_160209_OV3_P5_SK_JHa_373-Pe3_01.mzid_Mudpit_160209_OV3_P5_SK_JHa_373-Pe3_01.MGF
Mudpit_160213_OV3_P5_SK_JHa_372-Pe1_01.mzid_Mudpit_160213_OV3_P5_SK_JHa_372-Pe1_01.MGF
Mudpit_160215_OV3_P5_SK_JHa_372-Pe3_01.mzid_Mudpit_160215_OV3_P5_SK_JHa_372-Pe3_01.MGF
Mudpit_160214_OV3_P5_SK_JHa_372-Pe2_01.mzid_Mudpit_160214_OV3_P5_SK_JHa_372-Pe2_01.MGF
OP3_AA_20170330_fw_rev_cont.fasta
PXD025008_community_annotated.sdrf.tsv
files_overview_for_PRIDE_submission_20210325.pdf
"
for f in $FILES; do
  want=$(curl -sSIf -A "$UA" --max-time 60 "$BASE/$f" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r')
  for i in $(seq 1 30); do
    have=$(stat -f %z "$f" 2>/dev/null || echo 0)
    [ "$have" -lt 1000 ] && [ "$have" -gt 0 ] && rm -f "$f" && have=0      # an error page, not data
    [ -n "$want" ] && [ "$have" = "$want" ] && break
    curl -sSf -C - --connect-timeout 30 --max-time 3600 -A "$UA" -o "$f" "$BASE/$f" || sleep 10
  done
  echo "$f $(stat -f %z "$f" 2>/dev/null) of $want"
done
