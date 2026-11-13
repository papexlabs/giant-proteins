#!/usr/bin/env bash
# Build a clean release tree from the committed state (HEAD) and check it. Prepares only: it creates no
# git repository, no remote, and pushes nothing.
#
#   release/make_release.sh <target_dir> [--with-notes] [--verify]
#
#   --with-notes  include notes/ (working notes, first person; excluded by default, see CHECKLIST below)
#   --verify      in the release tree, restore the stored RNA-seq results and rerun the grading chain;
#                 the regenerated data/locus_checks.tsv, rnaseq_support.tsv and leaderboard.tsv must equal
#                 the committed ones (no network needed)
#
# What the tree contains: `git archive HEAD`, i.e. tracked files only. Untracked files (cloud/launched,
# cloud/stop_utc, cloud/watch.sh, cloud_out_watch.log), git-ignored files (cloud/aws.env, cloud/ip, cloud_out/,
# reads, raw spectra, genomes, tool binaries) are never included. The tree is then scanned for account-specific
# or personal strings; any hit stops the script with a non-zero exit.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[ $# -ge 1 ] || { echo "usage: $0 <target_dir> [--with-notes] [--verify]"; exit 1; }
DEST="$1"; shift
NOTES=0; VERIFY=0
for a in "$@"; do case "$a" in --with-notes) NOTES=1 ;; --verify) VERIFY=1 ;; *) echo "unknown option $a"; exit 1 ;; esac; done
[ -e "$DEST" ] && { echo "$DEST exists; choose a new directory"; exit 1; }
cd "$ROOT"
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "WARNING: uncommitted changes to tracked files are NOT in the release (it is built from HEAD)"
fi
mkdir -p "$DEST"
git archive --format=tar HEAD | tar -x -C "$DEST"
[ "$NOTES" = 1 ] || rm -rf "$DEST/notes"

echo "== scanning $DEST for account-specific or personal strings"
COMMON='/Users/|/home/[a-z]+/|kobon|sg-[0-9a-f]{8,17}|subnet-[0-9a-f]{8,17}|ami-[0-9a-f]{8,17}|vpc-[0-9a-f]{8,17}|\bi-[0-9a-f]{17}\b|aws_secret_access_key|BEGIN [A-Z ]*PRIVATE KEY|ghp_[A-Za-z0-9]{30,}|xox[abp]-[A-Za-z0-9-]+|[A-Za-z0-9._%+-]+@(gmail|outlook|hotmail|yahoo|icloud)\.[a-z]+|ahmet\.?tahir'
# /home/ubuntu is the generic user of the cloud image, not a person; this script holds the patterns themselves
HITS=$(cd "$DEST" && grep -rIoE "(${COMMON}|AKIA[0-9A-Z]{16})" . --exclude='*.fasta' --exclude='*.fa' \
       --exclude='make_release.sh' | grep -v ':/home/ubuntu/$' || true)
# protein FASTA files: AKIA-like runs of capitals occur in protein sequences, so only the other patterns apply
FAHITS=$(cd "$DEST" && grep -rIoE "(${COMMON})" . --include='*.fasta' --include='*.fa' | grep -v ':/home/ubuntu/$' || true)
IPS=$(cd "$DEST" && grep -rIohE '\b([0-9]{1,3}\.){3}[0-9]{1,3}\b' --include='*.sh' --include='*.py' --include='*.md' --include='*.txt' --include='*.json' --include='*.log' . \
      | grep -vE '^(0\.|127\.0\.0\.1|[0-9]\.[0-9]+\.[0-9]+\.[0-9]+$)' || true)
if [ -n "$HITS$FAHITS" ] || [ -n "$IPS" ]; then
  echo "FOUND (release tree not clean):"; echo "$HITS" | sort -u | head -30; echo "$FAHITS" | sort -u | head; echo "$IPS" | sort -u | head
  exit 2
fi
echo "clean: no home paths, AWS IDs, profile/key names, instance IDs, IPv4 addresses, keys, tokens or personal e-mail"

if [ "$VERIFY" = 1 ]; then
  echo "== verifying that the stored results reproduce the released grades"
  ( cd "$DEST" && mkdir -p cloud_out && tar xzf data/rnaseq_run_results.tar.gz -C cloud_out \
    && for f in locus_checks rnaseq_support leaderboard; do cp "data/$f.tsv" "/tmp/release_check_$f.$$"; done \
    && python3 scripts/13_grade.py > /dev/null && python3 cloud/summarize.py > /dev/null \
    && python3 scripts/15_apply_rnaseq.py > /dev/null && python3 scripts/20_apply_dna.py > /dev/null \
    && python3 scripts/16_leaderboard.py > /dev/null \
    && for f in locus_checks rnaseq_support leaderboard; do cmp -s "data/$f.tsv" "/tmp/release_check_$f.$$" \
         && echo "  data/$f.tsv reproduced" || { echo "  data/$f.tsv DIFFERS"; exit 3; }; rm -f "/tmp/release_check_$f.$$"; done \
    && rm -rf cloud_out )
fi

cat <<'EOF'
== CHECKLIST before any public release (prepared, not done)
[ ] Licence: choose one (release/LICENSE_PROPOSAL.md) and add LICENSE to the tree.
[ ] Git history: start a FRESH history in the release tree (git init; one initial commit). The working
    repository's history contains a personal e-mail address in every commit and, in older versions of
    cloud/launch.sh and cloud/sync.sh, the AWS profile, key-pair name, security-group, subnet and AMI IDs.
[ ] Commit identity for the public history: chosen by the author.
[ ] notes/: excluded by default (working notes in the first person, mention the coordinating agent and costs);
    include with --with-notes only if wanted. Three numbers in paper/NUMBERS.md cite notes/.
[ ] Statements in paper/manuscript.tex (author contributions, use of AI tools, competing interests, funding,
    acknowledgements) and the contact address: written by the author.
[ ] Zenodo: archive the tagged release; add the DOI to the manuscript's data availability section afterwards.
[ ] Third-party material: no tool binaries, no paper full texts, no reads or spectra are in the tree
    (see DATA_MANIFEST.tsv); proteomics/sycon/raw/Protein_Report_*.xls is the PXD060105 authors' deposited
    report (PRIDE), kept because analyze.py reads it.
EOF
echo "release tree: $DEST"
