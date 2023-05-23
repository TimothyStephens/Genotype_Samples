#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(

cat "${snakemake_input['samples']}" \
  | awk -F'\t' '{
    if(NR==1){
      SAMPLE_ID=0
      LIB_TYPE=0
      for(i=1; i<=NF; i++){
        if($i == "sample_id"){
          SAMPLE_ID=i
        }
        if($i == "lib_type"){
          LIB_TYPE=i
        }
      }
    }
    print $SAMPLE_ID"\t"$LIB_TYPE"\t"$0
  }' \
  | sort -k1,1 -k2,2 \
  | awk '!seen[$1]++{ if($1=="sample_id"){HEADER=$0}else{LINES[NR]=$0} } END {print HEADER; for(i in LINES){print LINES[i]}}' \
  | "${snakemake_params['add_values']}" \
    -c 1 -d "NA" \
    -a "${snakemake_input['groups']}" \
  | "${snakemake_params['add_values']}" \
    -c 1 -d "NA" \
    -a <(awk -F'\t' 'NR==1 || $2=="denoised"{print $1"\t"$13"\t"$14}' "${snakemake_input['nQuire']}") \
  | cut -f3- \
  > "${snakemake_output['samples']}"

) 1> "${snakemake_log[0]}" 2>&1


