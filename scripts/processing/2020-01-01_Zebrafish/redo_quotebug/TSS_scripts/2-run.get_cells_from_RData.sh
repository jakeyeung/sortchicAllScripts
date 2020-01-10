#!/bin/sh
# Jake Yeung
# 2-run.get_cells_from_RData.sh
#  
# 2019-11-08

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/Zebrafish/make_count_tables_TSS/get_cells_from_RData.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_ZFbonemarrow"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_PZ-ZF-All"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_PZ-ZF-All/good_cells"
[[ ! -d $outmain ]] && mkdir $outmain
outfmerged=$outmain/"good_cells_all_merged.csv"
[[ -e $outfmerged ]] && echo "$outfmerged found, exiting" && exit 1

for f in `ls -d $inmain/PZ-ZF*.rds`; do
  fbase=$(basename $f)
  fbase=${fbase%.*}
  outf=$outmain/$fbase.csv
  Rscript $rs $f $outf
  cat $outf >> $outfmerged
done

