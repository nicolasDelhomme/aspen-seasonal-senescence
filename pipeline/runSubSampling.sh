#!/bin/bash -l
#SBATCH -p core -n 1 -t 1:00:00 --mail-type=FAIL,END
set -eux

source ${SLURM_SUBMIT_DIR:-$(pwd)}/../UPSCb-common/src/bash/functions.sh

USAGETXT=\
"
Usage: $0 <input file> <number of reads> <output dir>

Note: The number of reads will be multipled by 4 to obtain the number of lines; i.e. we expect a gz fastq file as input

"

if [ $# -ne 3 ]; then
  abort "This script expects 3 parameters"
fi

if [ ! -f $1 ]; then
  abort "The input file does not exist"
fi

if [ $2 -lt 1 ]; then
  abort "You need to select at least one line"
fi

if [ ! -d $3 ]; then
  abort "The output directory does not exist"
fi

gunzip -c $1 | head -n $(expr $2 "*" 4) | gzip -c > $3/$(basename ${1/.fq.gz/_sub_$2.fq.gz/})
