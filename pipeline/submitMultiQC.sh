#!/bin/bash

## be verbose and print
set -ex

## source functions
source ../UPSCb-common/src/bash/functions.sh

## vars
proj=snic2018-8-303
mail=nicolas.delhomme@umu.se
in=/proj/uppoff2019006/uppstore2018168/P12869_results/
out=$in/multiqc

## tools
module load bioinfo-tools MultiQC

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## execute
sbatch -A $proj -t 2:00:00 --mail-user=$mail \
-e $out.err -o $out.out -J multiqc \
../UPSCb-common/pipeline/runMultiQC.sh $in $out

