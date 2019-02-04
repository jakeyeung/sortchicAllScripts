#!/bin/sh
# Jake Yeung
# 2019-02-04_send_bamlist.sh
# For merging bams 
# 2019-02-04

indir="/Users/yeung/projects/scchic/data/outputs/cell_clusters_bin_H3K4me1"
outdir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bamlists_for_merging/"

rsync -avrL --copy-links $indir $outdir
