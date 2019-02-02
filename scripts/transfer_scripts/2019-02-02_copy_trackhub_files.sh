#!/bin/sh
# Jake Yeung
# 2019-02-02_copy_trackhub_files.sh
# Copy trackhub try to visualize 
# 2019-02-02

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/trackhub_files/H3K4me1_peaks"
outdir="/Users/yeung/data/trackhubs/"

rsync -avrL --copy-links $indir $outdir
