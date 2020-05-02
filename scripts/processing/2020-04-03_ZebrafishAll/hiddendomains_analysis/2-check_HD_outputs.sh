#!/bin/sh
# Jake Yeung
# 2-check_HD_outputs.sh
# Check outputs have 25 chromosmoes in all 
# 2020-04-24
# 

outmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/hiddendomains_outputs"

for outdir in `ls -d ${outmain}/PZ*`; do
    bname=$(basename $outdir)
    hdout="$outdir/${bname}_analysis.bed"
    [[ ! -e $hdout ]] && echo "$hdout not found, exiting" && exit 1
    echo $bname
    cut -f1 $hdout  | sort -k1,1V | uniq -c | wc -l 
done
