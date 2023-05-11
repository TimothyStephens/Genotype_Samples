#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(
for FILE in "${snakemake_input[@]}";
do
  SAMPLE=$(basename ${FILE%.samtools_stats.txt})
  REF_NAME=$(basename $(dirname $(dirname $FILE)))
  awk -F'\t' \
      -vSAMPLE="${SAMPLE}" \
      -vREF_NAME="${REF_NAME}" \
  'BEGIN{
    NO_READS=0
    NO_MAPPED=0
  } $1=="SN" {
    if($2=="sequences:"){NO_READS=$3}
    if($2=="reads mapped:"){NO_MAPPED=$3}
  } END {
     printf "%s\t%s\t%.2f\n", SAMPLE, REF_NAME, (NO_MAPPED/NO_READS)*100
  }' "${FILE}"
done \
  | sort -k2,2 -k1,1 \
  | awk -F'\t' '{
    SAMPLE[$1]++
    REF[$2]++
    MATRIX[$1][$2]=$3
  } END {
    L="sample_id"
    for(R in REF){
      L=L"\t"R
    }
    print L
    
    for(S in SAMPLE){
      L=S
      for(R in REF){
        L=L"\t"MATRIX[S][R]
      }
      print L
    }
  }' > "${snakemake_output[0]}"
) 1> "${snakemake_log[0]}" 2>&1


