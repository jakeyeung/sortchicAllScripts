#!/bin/sh
# Jake Yeung
# calculate_mappability_summarize.sh
#  
# 2019-06-10

inmain="/hpc/hub_oudenaarden/jyeung/data"
cd $inmain

grep "total\|samtools\|mapped (" calculate_mappability_cat_across_3_experiments.out | sed "s/\/hpc\/hub_oudenaarden\/avo\/scChiC\/raw_demultiplexed\///g"

# # count total mapped counts
# totalmapped=`grep "total\|mapped (" calculate_mappability.out   | awk '(NR+1)%2==1' | cut -d" " -f1 | awk '{n += $1}; END{print n}'`
# 
# # count total counts 
# total=`grep "total\|mapped (" calculate_mappability.out   | awk 'NR%2==1' | cut -d" " -f1 | awk '{n += $1}; END{print n}'`
# 
# echo $totalmapped
# echo $total
# 
# awk -v var1=$totalmapped -v var2=$total 'BEGIN { print  ( var1 / var2 ) }'
# 
