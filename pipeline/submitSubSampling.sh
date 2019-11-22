#!/bin/bash -l

set -eux

proj=snic2019-8-310
mail=jenna.lihavainen@umu.se
MREADS=20000000

in=/proj/snic2019-30-28/senescence-RNA-Seq/results/trimmomatic
out=$in/subsample

if [ ! -d $out ]; then
  mkdir -p $out
fi

for f in $(find $in -name "P12255_1*_trimmomatic_[1,2].fq.gz"); 
do
  sbatch -A $proj -p core -n 1 -e $out/${f//.fq.gz/.err} \
  -o $out/${f//.fq.gz/.out} -t 1:00:00 runSubsampling.sh $f $MREADS $out
done
