#!/bin/sh
# Jake Yeung
# 2-rerun_scripts.sh
#  
# 2020-02-18

inf="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/4-run_motevo/log/scriptpaths_to_rerun.txt"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

while read row; do
    # echo "qsub $row"
    qsub $row
done < $inf
