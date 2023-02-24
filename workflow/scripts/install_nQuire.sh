#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

# https://github.com/clwgg/nQuire
(
cd "${CONDA_PREFIX}/lib"
rm -fr nQuire

git clone --recursive https://github.com/clwgg/nQuire
cd nQuire
make submodules
make
cp nQuire "${CONDA_PREFIX}/bin/"
) 1> "${snakemake_log[0]}" 2>&1
touch "${snakemake_output[0]}"


