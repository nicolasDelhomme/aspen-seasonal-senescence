#!/bin/bash -l

set -eux

proj=snic2017-7-326
mail=jenna.lihavainen@umu.se

in=/proj/uppoff2019006/senescence-RNA-Seq/results/raw
out=/proj/uppoff2019006/senescence-RNA-Seq/results
start=3
end=6

module load bioinfo-tools FastQC trimmomatic SortMeRNA/2.1b

if [ ! -d $out ]; then
  mkdir -p $out
fi


for f in $(find $in -name "*_1.fastq.gz"); 
do
  fnam=$(basename ${f/_1.fastq.gz/})
  bash ../UPSCb-common/pipeline/runRNASeqPreprocessing.sh -s $start -e $end \
  $proj $mail $f $in/${fnam}_2.fastq.gz $out
done
