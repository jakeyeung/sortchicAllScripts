#!/bin/sh
# Jake Yeung
# 2019-03-26_send_bamnames_to_server.sh
# Send it 
# 2019-03-26

indir="/Users/yeung/data/scchic/tables/bamlist_for_peak_analysis_build95"
outdir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook"

scp -r $indir $outdir
