#!/bin/bash
# EC2 user-data template for the giant-protein RNA-seq stage (runs once as root at first
# boot; Ubuntu 24.04 arm64, Graviton). launch.sh fills @@SHUTDOWN_MIN@@ and @@DEADLINE_EPOCH@@.
# 1. Arms the hard stop FIRST: the box powers off after @@SHUTDOWN_MIN@@ minutes whatever
#    happens; it is launched with instance-initiated shutdown = terminate, so power-off
#    terminates it and deletes its disk (DeleteOnTermination).
# 2. Writes the deadline for the driver's own stop rule (/etc/gp-deadline, epoch seconds).
# 3. Installs build tools and samtools (Ubuntu package), builds pinned minimap2 2.31 from
#    source (sha256-checked).
# 4. Writes /opt/gp/READY when done; launch.sh waits for it, then starts the driver.
# No credentials and no personal data on the box; downloads use a generic User-Agent and
# read only public data (NCBI genomes FTP over HTTPS, ENA FASTQ over HTTPS).
set -euxo pipefail
shutdown -h +@@SHUTDOWN_MIN@@ "giant-proteins hard stop"
echo @@DEADLINE_EPOCH@@ > /etc/gp-deadline
exec > /var/log/gp-user-data.log 2>&1
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update
apt-get -qq install -y build-essential zlib1g-dev python3 curl rsync samtools > /dev/null
UA='PapexResearchBot/1.0 (+https://papex.org)'
mkdir -p /opt/gp/src && cd /opt/gp/src
curl -fsSL -A "$UA" -o minimap2.tar.bz2 https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31.tar.bz2
echo "c1351de6319c123369c2f4f37ba0ccf18c7ace47e2b1c0a35e30056b4a3bd9c9  minimap2.tar.bz2" | sha256sum -c -
tar xjf minimap2.tar.bz2
(cd minimap2-2.31 && make -s -j"$(nproc)" arm_neon=1 aarch64=1 > /dev/null)
install -m 755 minimap2-2.31/minimap2 /usr/local/bin/minimap2
mkdir -p /home/ubuntu/gp /home/ubuntu/gpwork && chown -R ubuntu:ubuntu /home/ubuntu/gp /home/ubuntu/gpwork
echo "* soft nofile 65536" >> /etc/security/limits.conf
echo "* hard nofile 65536" >> /etc/security/limits.conf
samtools --version | head -1 > /opt/gp/samtools.version
minimap2 --version > /opt/gp/READY
