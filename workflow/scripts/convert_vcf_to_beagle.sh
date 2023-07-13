#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(
GZVCF="${snakemake_input['vcf']}"
BEAGLE="${snakemake_output['beagle']}"

while read SEQID;
do
  vcftools --gzvcf "${GZVCF}" --stdout --BEAGLE-PL --chr "${SEQID}"
done < <(zcat "${GZVCF}" | awk '$1~"^##contig="' | sed -e 's/.*ID=\([^,]*\).*/\1/') \
  | awk -F'\t' 'NR==1 || $1!="marker"' \
  | gzip -c \
  > "${BEAGLE}"

) 1> "${snakemake_log[0]}" 2>&1


