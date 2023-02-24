#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

# http://www.popgen.dk/software/index.php/PCAngsd
(
cd "${CONDA_PREFIX}/lib"
rm -fr pcangsd

git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd
python setup.py build_ext --inplace
pip3 install -e .
) 1> "${snakemake_log[0]}" 2>&1
touch "${snakemake_output[0]}"


