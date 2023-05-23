#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate Genotype_Samples; set -eu



## Envs
TAG="0.0.1"

USER="timothystephens"
PROG="genotype_samples"
DOCKERFILE="docker/Dockerfile"
BUILD_LOG="docker/build.log"



## Snakemake containerize
cd ../
./Genotype_Samples.py --module genotyping --configfile tests/config/config.Genotyping_small.yaml --containerize \
  | sed -e '/Running snakemake command/d' -e '/module=all/d' \
  > $DOCKERFILE



## Build docker image
docker build --no-cache -t $PROG:$TAG $DOCKERFILE 2>&1 | tee $BUILD_LOG



## Build docker image and push to repo
BUILD=$(awk '$1=="Successfully" && $2=="built" {print $3}' $BUILD_LOG)
echo "BUILD=${BUILD}"

docker tag $BUILD $USER/$PROG:${TAG}
docker tag $BUILD $USER/$PROG:${TAG%.*}
docker tag $BUILD $USER/$PROG:${TAG%%.*}
docker tag $BUILD $USER/$PROG:latest

docker push --all-tags $USER/$PROG



## Cleanup docker and test files
docker image prune -af
docker image rm -f $(docker image ls | awk 'NR>1{print $3}')
# rm -fr .cache .keras .parallel .snakemake resources results


