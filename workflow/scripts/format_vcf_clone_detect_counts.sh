#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(
awk -F',' 'FNR==NR{
  names[FNR]=$1
} FNR!=NR {
  if(FNR==1){
    for(i=1; i<=length(names); i++){
      for(j=1; j<=length(names); j++){
        if(names[i] == names[j]){
          M[names[i]][names[j]]="100"
        } else {
          M[names[i]][names[j]]="NA"
        }
      }
    }
  } else { 
    M[$1][$2]=$7
    M[$2][$1]=$7
  }
} END {
  L="sample_id"; 
  for(i=1; i<=length(names); i++){
    L=L"\t"names[i]
  } 
  print L
  for(i=1; i<=length(names); i++){
    L=names[i]
    for(j=1; j<=length(names); j++){
      L=L"\t"M[names[i]][names[j]]
    }
    print L
  } 
}' "${snakemake_input[0]}" "${snakemake_input[1]}" > "${snakemake_output[0]}"
) 1> "${snakemake_log[0]}" 2>&1


