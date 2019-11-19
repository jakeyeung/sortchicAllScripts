#!/bin/sh
# Jake Yeung
# 3-run.RData_to_rds_for_inmat.sh
#  
# 2019-09-30

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells"
outdir=$indir

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_stem_cells/1-run_LDA/RData_to_rds_for_inmat.R"

for inf in `ls -d $indir/*.RData`; do
    bbase=$(basename $inf)
    bbase=${bbase%.*}
    echo $bbase
    outf=$outdir/$bbase.rds
    Rscript $rs $inf $outf
done

