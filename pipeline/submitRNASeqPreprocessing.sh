#!/bin/bash -l

set -eux

proj=snic2018-8-303
mail=nicolas.delhomme@umu.se

in=/proj/uppstore2018168/P12869
out=/proj/uppoff2019006/uppstore2018168/P12869_results
start=3
end=6

module load bioinfo-tools FastQC trimmomatic SortMeRNA/2.1b

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in $(find $in -name "*_1.fastq.gz"); 
do
  fnam=$(basename ${f/_1.fastq.gz/})
  bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end \
  $proj $mail $f $in/${fnam}_2.fastq.gz $out
done
