#!/bin/sh
# Jake Yeung
# run.run.load_mara_save_output_array.sh
#  
# 2019-04-10

# goes from 1 to 23170

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/run.load_mara_save_output_array.sh"

[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

qsub -t 1-23170 $bs
