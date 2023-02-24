#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

# https://github.com/clwgg/nQuire
(
awk 'BEGIN{OFS=FS="\t"} {
  if(NR==1){
    $1="sample_id\ttype";
    print
  } else if($1~".denoised.bin$"){
    split($1,a,"/");
    b=a[length(a)];
    sub(".denoised.bin$","",b);
    $1=b"\tdenoised";
    print;
  } else if($1~".bin$"){
    split($1,a,"/");
    b=a[length(a)];
    sub(".bin$","",b);
    $1=b"\tnormal";
    print;
  } 
}' ${snakemake_input[0]} > ${snakemake_output[0]}
) 1> "${snakemake_log[0]}" 2>&1


