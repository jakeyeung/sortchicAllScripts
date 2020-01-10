#!/bin/sh
# Jake Yeung
# 7b-link_retagged_bams_together.sh
#  
# 2019-12-01

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all/raw_demultiplexed.first_try"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag"

for inf in `ls -d $inmain/PZ*/tagged.retagged/*bam*`; do
    echo $inf
    infbase=$(basename $inf)
    outf=$outdir/$inbase
    ln -s $inf $outf
done
