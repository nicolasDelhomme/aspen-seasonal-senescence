#!/bin/bash -l
out=$(realpath ../data/seidr/aggregated)
mail=kristina.benevides@umu.se

module load bioinfo-tools seidr-devel

sbatch -o $out/aggregated.out -e $out/aggregated.err --mail-user $mail ../UPSCb-common/pipeline/runSeidrAggregate.sh $out $out/*.sf

