#!/bin/bash -l

## be verbose and print
set -ex

proj=snic2017-7-326
mail=jenna.lihavainen@umu.se

## process the argument (just change in and out)
in=/proj/uppoff2019006/senescence-RNA-Seq/results/trimmomatic
ref=/proj/uppoff2019006/indices/Potra02/Potra02_v2dot2_transcripts_salmon-v14dot1.inx
out=/proj/uppoff2019006/senescence-RNA-Seq/results/salmon
bind=/proj/uppoff2019006:/proj/uppoff2019006
img=/proj/uppoff2019006/singularity/singularity/salmon-0.14.1.simg

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## for every file
for f in $(find $in -name "*sortmerna_trimmomatic_1.fq.gz"); do
  fnam=$(basename ${f/_1.fq.gz/})

  ## execute
  sbatch -A $proj --mail-user=$mail \
  -e $out/$fnam.err -o $out/$fnam.out -J salmon.$fnam \
  ../UPSCb-common/pipeline/runSalmon.sh -b $bind \
  -i $img $ref $f $in/${fnam}_2.fq.gz $out

done
