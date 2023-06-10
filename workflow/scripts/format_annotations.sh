#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(

# Cleanup samples.tsv file. Remove commented and blank lines
awk -F'\t' '$1!~"^#" && $0!=""' \
    "${snakemake_input['samples']}" \
  > "${snakemake_output['samples']}"

# Get indexes of columns that we want to use for annotations in plots
COL_INDEX=$(
  awk -F'\t' '{
    if(NR==1){
      for(i=1; i<=NF; i++){
        if($i != "sample_id" && $i != "unit" && $i != "fq1" && $i != "fq2"){
          H[i]++
        }
      }
    }
  }END{
    for(i in H){
      print i
    }
  }' "${snakemake_output['samples']}"
)

# Add colors to all unique values in the annotation columns
awk -F'\t' '{ if(FNR==NR){ COL[FNR]=$1 }else{ print $1"\t"COL[FNR] } }' \
    "${snakemake_input['color_list']}" \
    <(
      while read idx;
      do
        awk -F'\t' -vidx="$idx" 'NR>1{print $idx}' "${snakemake_output['samples']}" | sort | uniq
      done < <(echo -e "$COL_INDEX") \
    ) \
  > "${snakemake_output['color_list']}"

# Append ploidy colors to list
echo -e "Diploid\t#1b9e77\nTriploid\t#d95f02\nTetraploid\t#7570b3\nUnknown\t#808080" \
  >> "${snakemake_output['color_list']}"

# Append "Group" IDs to end of list
awk -F'\t' '{ print "Group"NR"\t"$1 }' \
    "${snakemake_input['color_list']}" \
  >> "${snakemake_output['color_list']}"
echo -e "Ungrouped\t#808080" \
  >> "${snakemake_output['color_list']}"

) 1> "${snakemake_log[0]}" 2>&1


