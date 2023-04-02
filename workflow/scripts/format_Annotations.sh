#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(
"${snakemake_params['add_values']}" \
  -c 2 \
  -i "${snakemake_input['groups']}" \
  -a <(echo -e "group_id\tgroup_id_color\nGroup1\t#a6cee3\nGroup2\t#1f78b4\nGroup3\t#b2df8a\nGroup4\t#33a02c\nGroup5\t#fb9a99\nGroup6\t#e31a1c\nGroup7\t#fdbf6f\nGroup8\t#ff7f00\nGroup9\t#cab2d6\nGroup10\t#6a3d9a\nGroup11\t#ffff99\nGroup12\t#b15928\nGroup13\t#8dd3c7\nGroup14\t#ffffb3\nGroup15\t#bebada\nGroup16\t#fb8072\nGroup17\t#80b1d3\nGroup18\t#fdb462\nGroup19\t#b3de69\nGroup20\t#fccde5\nGroup21\t#d9d9d9\nGroup22\t#bc80bd\nGroup23\t#ccebc5\nGroup24\t#ffed6f\nUngroup\t#808080") \
  -d "#E5E4E2" \
| "${snakemake_params['add_values']}" \
  -a <(\
    "${snakemake_params['add_values']}" \
        -c 2 \
        -i <(awk -F'\t' 'NR==1 || $2=="denoised"{print $1"\t"$13}' "${snakemake_input['nQuire']}") \
        -a <(echo -e "best_ploidy_model\tbest_ploidy_model_color\nDiploid\t#1b9e77\nTriploid\t#d95f02\nTetraploid\t#7570b3")
        -d "#E5E4E2" \
  ) \
  > "${snakemake_output[0]}"
) 1> "${snakemake_log[0]}" 2>&1


