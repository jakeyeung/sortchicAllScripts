#!/bin/sh
# Jake Yeung
# 2019-03-29_send_bed_files.sh
# Send bed files
# 2019-03-29

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output_build95"
outdir="/Users/yeung/data/scchic/from_cluster/hiddenDomains_output"

scp -r $indir $outdir
