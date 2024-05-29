#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(
awk -F' ' 'FNR==NR{
  names[FNR]=$1
} FNR!=NR {
  if(FNR==1){
    L="sample_id"; 
    for(i=1; i<=NF; i++){
      L=L"\tP"i
    };
    print L
  };
  gsub(" ","\t",$0)
  print names[FNR]"\t"$0
}' "${snakemake_input[0]}" "${snakemake_input[1]}" > "${snakemake_output[0]}"
) 1> "${snakemake_log[0]}" 2>&1


