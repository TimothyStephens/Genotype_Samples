#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

(
cp "${snakemake_input[0]}" "${snakemake_output[0]}"
) 1> "${snakemake_log[0]}" 2>&1


