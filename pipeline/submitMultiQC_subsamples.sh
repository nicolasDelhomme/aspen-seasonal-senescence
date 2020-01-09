#!/bin/bash

## be verbose and print
set -ex

## source functions
source ../UPSCb-common/src/bash/functions.sh

## vars
proj=snic2019-8-324
mail=jenna.lihavainen@umu.se
in=/proj/snic2019-30-28/senescence-RNA-Seq/results/
out=$in/multiqc/subsample

## tools
module load bioinfo-tools MultiQC

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## execute
sbatch -n 4 -A $proj -t 2:00:00 --mail-user=$mail \
-e $out.err -o $out.out -J multiqc \
../UPSCb-common/pipeline/runMultiQC.sh $in $out

