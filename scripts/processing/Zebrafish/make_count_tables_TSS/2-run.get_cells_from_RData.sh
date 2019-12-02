#!/bin/sh
# Jake Yeung
# 2-run.get_cells_from_RData.sh
#  
# 2019-11-08

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/Zebrafish/make_count_tables_TSS/get_cells_from_RData.R"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_ZFbonemarrow"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_ZFbonemarrow/good_cells"
outfmerged=$outmain/"good_cells_all_merged.csv"

[[ -e $outfmerged ]] && echo "$outfmerged found, exiting" && exit 1

for f in `ls -d $inmain/ZF-H3*.2019-11-04.RData`; do
  fbase=$(basename $f)
  fbase=${fbase%.*}
  outf=$outmain/$fbase.csv
  Rscript $rs $f $outf
done

cat $outf >> $outfmerged
