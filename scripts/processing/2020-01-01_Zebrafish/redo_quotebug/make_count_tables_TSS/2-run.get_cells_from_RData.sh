#!/bin/sh
# Jake Yeung
# 2-run.get_cells_from_RData.sh
#  
# 2019-11-08

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/Zebrafish/make_count_tables_TSS/get_cells_from_RData.R"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_ZF"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_ZF/good_cells"
[[ ! -d $outmain ]] && mkdir $outmain
outfmerged=$outmain/"good_cells_all_merged.csv"
[[ -e $outfmerged ]] && echo "$outfmerged found, exiting" && exit 1

for f in `ls -d $inmain/ZF_AllMerged_*2019-12-05.rds`; do
  fbase=$(basename $f)
  fbase=${fbase%.*}
  outf=$outmain/$fbase.csv
  . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $f $outf
  cat $outf >> $outfmerged
done

