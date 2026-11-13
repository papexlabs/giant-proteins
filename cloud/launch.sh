#!/bin/bash
# Launch ONE spot box for the giant-protein RNA-seq stage, push code + job bundle, start
# the driver. eu-north-1 only (the only region the account may use).
#
#   usage: cloud/launch.sh <HARD_STOP_UTC>        e.g. cloud/launch.sh 2026-11-14T06:00Z
#          STAGE=followup cloud/launch.sh <HARD_STOP_UTC>
#            follow-up stage: DNA re-mapping checks (bundle_dna, dna_driver.py) then the RNA
#            re-run of the three collision loci (bundle_followup_rna), via followup_box.sh;
#            a driver.log that already ends an earlier stage ("powering off") is rotated first
#          STAGE=rep cloud/launch.sh <HARD_STOP_UTC>
#            RNA-seq of the repeat-flag models (bundle_rep_rna, 13 Nov): ONLY the RNA driver runs,
#            with --download (each FASTQ file to disk with resumed retries, deleted after its
#            run's result is written) and no saturation stop; the driver powers the box off when
#            its queue is done (rules fixed in cloud/bundle_rep_rna/RULES.txt before launch)
#
# Safety:
# - at most 32 vCPU (r8g.8xlarge first; m8g/r7g/c8g/c7g.8xlarge fallbacks), one-time spot,
#   MaxPrice per candidate, instance-initiated shutdown = terminate, root gp3 400 GB (500 MB/s) with
#   DeleteOnTermination, tag Name=giant-proteins-rna on instance and volume;
# - refuses to launch if a giant-proteins-rna box already runs, if running spot vCPU + 32
#   would leave fewer than 12 of the 300-vCPU spot quota free, or if the hard stop is
#   < 2 h or > 72 h away; the hard stop is armed on the box (shutdown -h +MIN) before
#   anything else, and the driver stops starting species it cannot finish in time;
# - never touches other instances (no stop/terminate/modify calls in this script).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
# account-specific settings (profile, key pair, security group, subnets, AMI) live in cloud/aws.env (git-ignored)
[ -f "$HERE/aws.env" ] || { echo "cloud/aws.env missing: copy cloud/aws.env.example and fill it in"; exit 1; }
set -a; . "$HERE/aws.env"; set +a
export AWS_PAGER=""
cd "$HERE"
[ $# -eq 1 ] || { echo "usage: $0 <HARD_STOP_UTC e.g. 2026-11-14T06:00Z>"; exit 1; }
STOP="$1"
read -r MIN EPOCH < <(python3 - "$STOP" <<'PY'
import sys, datetime as d
t = d.datetime.strptime(sys.argv[1].replace("Z", "+0000"), "%Y-%m-%dT%H:%M%z")
now = d.datetime.now(d.timezone.utc)
print(int((t - now).total_seconds() // 60), int(t.timestamp()))
PY
)
[ "$MIN" -ge 120 ] || { echo "hard stop is less than 2 h away ($MIN min): not launching"; exit 1; }
[ "$MIN" -le 4320 ] || { echo "hard stop is more than 72 h away ($MIN min): not launching (budget rule)"; exit 1; }
STAGE="${STAGE:-rna}"
[ -f bundle/jobs.tsv ] || { echo "cloud/bundle/jobs.tsv missing: run cloud/make_jobs.py first"; exit 1; }
case "$STAGE" in rna|followup|rep) ;; *) echo "unknown STAGE=$STAGE (rna, followup or rep)"; exit 1 ;; esac
if [ "$STAGE" = rep ]; then
  [ -f bundle_rep_rna/jobs.tsv ] && [ -f bundle_rep_rna/RULES.txt ] \
    || { echo "bundle_rep_rna/jobs.tsv or RULES.txt missing"; exit 1; }
fi
if [ "$STAGE" != rna ]; then
  [ "$STAGE" = rep ] || { ls bundle_dna/jobs/*.json > /dev/null && [ -f bundle_followup_rna/jobs.tsv ]; } \
    || { echo "follow-up bundles missing (bundle_dna/jobs, bundle_followup_rna)"; exit 1; }
  if grep -q "powering off" "$ROOT/cloud_out/driver.log" 2>/dev/null; then
    mv "$ROOT/cloud_out/driver.log" "$ROOT/cloud_out/driver.$(date -u +%Y%m%dT%H%M%SZ).log"
    echo "rotated the finished stage's driver.log"
  fi
fi

# guard 1: no second box of this stage
EXIST=$(aws ec2 describe-instances --filters Name=tag:Name,Values=giant-proteins-rna \
  Name=instance-state-name,Values=pending,running --query 'Reservations[].Instances[].InstanceId' --output text)
[ -z "$EXIST" ] || { echo "a giant-proteins-rna box is already running: $EXIST"; exit 1; }
# guard 2: keep >= 12 spot vCPU of the 300 quota free (other project's boxes are counted, not touched)
USED=$(aws ec2 describe-instances --filters Name=instance-lifecycle,Values=spot \
  Name=instance-state-name,Values=pending,running \
  --query 'Reservations[].Instances[].[CpuOptions.CoreCount,CpuOptions.ThreadsPerCore]' --output text \
  | awk '{s+=$1*$2} END{print s+0}')
[ $((USED + 32)) -le 288 ] || { echo "spot vCPU in use $USED: launching 32 would leave < 12 free"; exit 1; }
echo "spot vCPU in use: $USED; after launch: $((USED + 32)) of 300"

sed "s/@@SHUTDOWN_MIN@@/$MIN/; s/@@DEADLINE_EPOCH@@/$EPOCH/" bootstrap.sh > user-data.rendered.sh
A=$GP_SUBNET_A; B=$GP_SUBNET_B; C=$GP_SUBNET_C   # three availability zones (cloud/aws.env)
# candidates "type:subnet:maxprice" (USD/h), ordered by what was least likely to be reclaimed
# (read-only checks on 12 Nov 2026, see notes section 13): AWS spot placement score (1-10)
# for one instance per AZ was 3 for r8g/m8g/c8g/c7g.8xlarge in every AZ and for r7g.8xlarge
# only in eu-north-1b (2 in 1a/1c; the first run was reclaimed 4 times in 7 h on r7g in 1c);
# Spot Advisor (region level) puts all of them in the <5 % interruption band. 7-day spot
# prices: r8g 0.20-0.26, m8g 0.18-0.30, r7g(1b) 0.18-0.20, c8g 0.17-0.25, c7g 0.15-0.23.
# r8g first (Graviton4, 256 GB: index builds of the largest genomes); m7g is left out
# because other workloads in the same account used m7g capacity at the time.
CANDS="r8g.8xlarge:$B:0.30 r8g.8xlarge:$A:0.30 r8g.8xlarge:$C:0.30 m8g.8xlarge:$A:0.32 r7g.8xlarge:$B:0.25 m8g.8xlarge:$B:0.32 c8g.8xlarge:$B:0.27 c7g.8xlarge:$B:0.25"
IID=""
for c in $CANDS; do
  T=${c%%:*}; r=${c#*:}; SUB=${r%%:*}; MP=${r#*:}
  OUT=$(aws ec2 run-instances --image-id "$GP_AMI" --instance-type "$T" --count 1 \
    --key-name "$GP_KEY_NAME" --security-group-ids "$GP_SECURITY_GROUP" --subnet-id "$SUB" \
    --instance-market-options "MarketType=spot,SpotOptions={MaxPrice=$MP,SpotInstanceType=one-time,InstanceInterruptionBehavior=terminate}" \
    --instance-initiated-shutdown-behavior terminate \
    --block-device-mappings '[{"DeviceName":"/dev/sda1","Ebs":{"VolumeSize":400,"VolumeType":"gp3","Throughput":500,"DeleteOnTermination":true}}]' \
    --user-data file://user-data.rendered.sh \
    --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=giant-proteins-rna}]' \
                         'ResourceType=volume,Tags=[{Key=Name,Value=giant-proteins-rna}]' \
    --query 'Instances[0].InstanceId' --output text 2>&1) || true
  case "$OUT" in i-*) IID=$OUT; echo "launched $T in $SUB, max \$$MP/h: $IID"; break ;;
    *) echo "no launch $T $SUB: $(echo "$OUT" | cut -c1-160)" ;; esac
done
[ -n "$IID" ] || { echo "LAUNCH FAILED (no capacity at the max prices)"; exit 1; }
aws ec2 wait instance-running --instance-ids "$IID"
IP=$(aws ec2 describe-instances --instance-ids "$IID" --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
echo "$IID $(date -u +%FT%TZ) $IP $T $MP stop=$STOP" >> launched; echo "$IP $IID" > ip
SSH="ssh -i $GP_KEY_FILE -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=20"
for i in $(seq 1 60); do $SSH -n ubuntu@"$IP" 'test -f /opt/gp/READY' 2>/dev/null && break; sleep 15; done
$SSH -n ubuntu@"$IP" 'test -f /opt/gp/READY' || { echo "bootstrap not ready after 15 min; box keeps its hard stop"; exit 1; }
rsync -az -e "$SSH" "$HERE/driver.py" "$HERE/rnaseq_locus.py" "$HERE/sam2fq.py" "$HERE/prep_genome.py" "$HERE/bundle" ubuntu@"$IP":/home/ubuntu/gp/
if [ "$STAGE" = followup ]; then
  rsync -az -e "$SSH" "$HERE/dna_check.py" "$HERE/dna_driver.py" "$HERE/followup_box.sh" "$HERE/bundle_dna" \
    "$HERE/bundle_followup_rna" ubuntu@"$IP":/home/ubuntu/gp/
elif [ "$STAGE" = rep ]; then
  rsync -az -e "$SSH" "$HERE/bundle_rep_rna" ubuntu@"$IP":/home/ubuntu/gp/
fi
# resume: results synced from an earlier (e.g. spot-interrupted) box are pushed back first
[ -d "$ROOT/cloud_out/results" ] && rsync -az -e "$SSH" "$ROOT/cloud_out/results/" ubuntu@"$IP":/home/ubuntu/gpwork/results/
# the synced driver.log goes back too, so the new box appends to it instead of starting over
[ -f "$ROOT/cloud_out/driver.log" ] && rsync -az -e "$SSH" "$ROOT/cloud_out/driver.log" ubuntu@"$IP":/home/ubuntu/gpwork/driver.log
if [ "$STAGE" = followup ]; then
  START="bash /home/ubuntu/gp/followup_box.sh"
elif [ "$STAGE" = rep ]; then
  START="python3 /home/ubuntu/gp/driver.py /home/ubuntu/gp/bundle_rep_rna /home/ubuntu/gpwork --threads 32 --par 8 --sat-cov 0 --download"
else
  START="python3 /home/ubuntu/gp/driver.py /home/ubuntu/gp/bundle /home/ubuntu/gpwork --threads 32 --par 8"
fi
$SSH -n ubuntu@"$IP" "mkdir -p /home/ubuntu/gpwork; (nohup setsid $START >> /home/ubuntu/gpwork/driver.out 2>&1 < /dev/null &); \
  sleep 3; pgrep -f 'driver.py|followup_box' | head -1; sudo shutdown --show 2>&1 | head -1"
echo "STARTED $STAGE stage on $IP ($IID). Sync results: cloud/sync.sh   Hard stop: $STOP"
