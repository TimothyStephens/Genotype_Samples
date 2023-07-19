#!/usr/bin/env bash

## Cant use pipefail becuse vcftools will return errorcode > 0 if no variants identified for a scaffold.
## This happens a lot with short or weird scaffolds, so is an expected behavior.
## vcftools error:
##	Error: Require GL or PL FORMAT tags in VCF file to output BEAGLE input.
#set -euo pipefail
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


