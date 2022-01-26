#!/bin/bash -l

account="u2019019"

# Load the tools
module load bioinfo-tools seidr-devel

# process the argument
aggregate=$(realpath ../data/seidr/aggregated/aggregated.sf)
backbone=$(realpath ../data/seidr/backbone)
out=$(realpath ../data/seidr/roc)
gs_potra02_pos=$(realpath ../data/seidr/gold_standard/Potra02_KEGG-based-positive-gold-standard.tsv)
gs_potra02_neg=$(realpath ../data/seidr/gold_standard/Potra02_KEGG-based-negative-gold-standard.tsv)

if [ ! -d $out ]; then
  mkdir -p $out
fi

# submit
for j in {1..10}; do
# for backbone files
  sbatch -A $account -o $out/${j}-percent.out \
  -e $out/${j}-percent.err -J roc-${j}  \
  ./runSeidrRoc.sh $backbone/backbone-${j}-percent.sf \
 $gs_potra02_pos $gs_potra02_neg $out/${j}-percent.roc
done

# for aggregated files
  sbatch -A $account -o $out/aggregated.out \
  -e $out/aggregated.err -J roc-aggr  \
  ./runSeidrRoc.sh $aggregate $gs_potra02_pos \
 $gs_potra02_neg $out/aggregated.roc
