#!/bin/bash
# Pull the small per-run results and the driver log from the box to the Mac
# (cloud_out/). Read-only on AWS (one describe call to check the box state).
#   usage: cloud/sync.sh            (uses cloud/ip written by launch.sh)
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
set -a; . "$HERE/aws.env"; set +a          # account-specific settings (git-ignored)
export AWS_PAGER=""
read -r IP IID < "$HERE/ip"
STATE=$(aws ec2 describe-instances --instance-ids "$IID" --query 'Reservations[0].Instances[0].State.Name' --output text 2>/dev/null)
echo "box $IID ($IP): $STATE"
SSH="ssh -i $GP_KEY_FILE -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=15"
mkdir -p "$ROOT/cloud_out"
rsync -az -e "$SSH" ubuntu@"$IP":/home/ubuntu/gpwork/results/ "$ROOT/cloud_out/results/" && \
rsync -az -e "$SSH" ubuntu@"$IP":/home/ubuntu/gpwork/driver.log "$ROOT/cloud_out/" && \
echo "synced: $(find "$ROOT/cloud_out/results" -name DONE | wc -l | tr -d ' ') species done, \
$(find "$ROOT/cloud_out/results" -name '*.json' | wc -l | tr -d ' ') run results" && \
tail -3 "$ROOT/cloud_out/driver.log"
