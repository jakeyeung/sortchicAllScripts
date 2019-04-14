#!/bin/sh
# Jake Yeung
# 3-npz_to_textfile.sh
#  
# 2019-04-01

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/public_data/lib/npz_to_textfile.py"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/Multimark_binsize-100000_comparison_build95"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/Multimark_binsize-100000_comparison_build95"

n=0
maxjobs=4
for inf in `ls -d $inmain/*.npz`; do
    inbase=$(basename $inf)
    inbase=${inbase%%.*}
    outf=$outmain/$inbase.txt
    python $ps $inf $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
