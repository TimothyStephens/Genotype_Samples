#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

# https://github.com/ANGSD/NgsRelate
(
export R_LIBS="${CONDA_PREFIX}/lib/R/library"

cd "${CONDA_PREFIX}/lib"
rm -fr genomescope2.0 KMC

URL="https://github.com/tbenavi1/genomescope2.0.git"
DIR="genomescope2.0"
git clone "$URL"
cd "$DIR/"
Rscript install.R
cd "${CONDA_PREFIX}/bin/"
ln -s "../lib/$DIR/genomescope.R"

cd ../

URL="https://github.com/tbenavi1/KMC.git"
DIR="KMC"
git clone "$URL"
cd "$DIR/"
make
cp bin/* "${CONDA_PREFIX}/bin/"
) 1> "${snakemake_log[0]}" 2>&1
touch "${snakemake_output[0]}"


