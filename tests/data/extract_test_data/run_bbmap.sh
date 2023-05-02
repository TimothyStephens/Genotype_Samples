#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="/home/timothy/programs/bbmap:$PATH"

REF="target_region.fa"
NCPUS=90

#### Start Script
while read LINE;
do
  ID=$(echo -e "$LINE" | awk -F'\t' '{print $1}')
  R1=$(echo -e "$LINE" | awk -F'\t' '{print $2}')
  R2=$(echo -e "$LINE" | awk -F'\t' '{print $3}')
  run_cmd "bbmap.sh threads=${NCPUS} in=${R1} in2=${R2} out=stdout.sam nodisk | ./filter_reads_for_bin_reassembly.py ../data ."
  run_cmd "mv ref.R1.fastq.gz ../data/${ID}_R1.fq.gz"
  run_cmd "mv ref.R2.fastq.gz ../data/${ID}_R2.fq.gz"
done < reads.tsv



