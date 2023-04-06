#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

# https://github.com/clwgg/nQuire
(
awk 'BEGIN{OFS=FS="\t"} {
  if(NR==1){
    $1="sample_id\ttype";
    print $0"\tbest_ploidy_model\tbest_ploidy_model_num"
  } else {
    if     ($6<$7 && $6<$8){BEST="Diploid\t2"}
    else if($7<$6 && $7<$8){BEST="Triploid\t3"}
    else if($8<$6 && $8<$7){BEST="Tetraploid\t4"}
    else {BEST="Unknown\tNA"}
    if($1~".denoised.bin$"){
      split($1,a,"/");
      b=a[length(a)];
      sub(".denoised.bin$","",b);
      $1=b"\tdenoised";
      print $0"\t"BEST;
    } else {
      split($1,a,"/");
      b=a[length(a)];
      sub(".bin$","",b);
      $1=b"\tnormal";
      print $0"\t"BEST;
    }
  } 
}' "${snakemake_input[0]}" \
  | "${snakemake_params['add_values']}" \
      -a "${snakemake_input[1]}" \
  | awk -F'\t' '{print $1"\t"$2"\t"$12"\t"$13"\t"$14"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' \
  > "${snakemake_output[0]}"
) 1> "${snakemake_log[0]}" 2>&1


