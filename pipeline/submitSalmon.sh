#!/bin/bash -l

## be verbose and print
set -ex

proj=snic2018-8-303
mail=nicolas.delhomme@umu.se

## process the argument
in=proj/uppoff2019006/uppstore2018168/P12869_results/trimmomatic
ref=/proj/uppoff2019006/indices/Potra01/Potra01-mRNA_salmon-v14dot1.idx
out=proj/uppoff2019006/uppstore2018168/P12869_results/Salmon
bind=/proj/uppoff2019006:/proj/uppoff2019006
img=//proj/uppoff2019006/singularity/singularity/salmon-0.14.1.simg

## check vars
if [ -z $UPSCb ]; then
    abort "The UPSCb var needs to be set."
fi

## load the tool - we actually use singularity
#module load bioinfo-tools salmon

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
  $UPSCb/pipeline/runSalmon.sh -b $bind \
  -i $img $ref $f $in/${fnam}_2.fq.gz $out

done
