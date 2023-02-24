#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

# https://github.com/ANGSD/NgsRelate
(
cd "${CONDA_PREFIX}/lib"
rm -fr htslib ngsRelate

git clone --recursive https://github.com/SAMtools/htslib
cd htslib/
make

cd ../

git clone https://github.com/ANGSD/ngsRelate
cd ngsRelate
make HTSSRC=../htslib/
cp ngsRelate "${CONDA_PREFIX}/bin/"
) 1> "${snakemake_log[0]}" 2>&1
touch "${snakemake_output[0]}"


