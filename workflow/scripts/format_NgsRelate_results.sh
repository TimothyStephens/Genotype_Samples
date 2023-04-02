#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(
cat "${snakemake_input[0]}" \
  | cut -f3- \
  | sed -e '1 s/ida/sample_id_1/' \
        -e '1 s/idb/sample_id_b/' \
  > "${snakemake_output[0]}"
) 1> "${snakemake_log[0]}" 2>&1


