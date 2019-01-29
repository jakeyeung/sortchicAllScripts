#!/bin/sh
# Jake Yeung
# 0-make_wms_list.sh
# Make WMs list for chekinv all motifs 
# 2019-01-29

inmain="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm10_weight_matrices_v2_split"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm10_v2_WMs.list"

cd $inmain

[[ -e $outf ]] && echo "$outf found, exiting for safety" && exit 1

ls -d *.pwm >> $outf
