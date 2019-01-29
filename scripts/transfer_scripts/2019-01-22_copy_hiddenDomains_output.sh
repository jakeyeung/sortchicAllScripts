#!/bin/sh
# Jake Yeung
# 2019-01-22_copy_hiddenDomains_output.sh
# Copy HD output 
# 2019-01-22

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output"
outdir="/tmp/hiddenDomains_out"

scp -r t2:$indir $outdir
