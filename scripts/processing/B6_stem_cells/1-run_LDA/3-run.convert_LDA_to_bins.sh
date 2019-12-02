#!/bin/sh
# Jake Yeung
# 3-run.convert_LDA_to_bins.sh
#  
# 2019-09-30

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_stem_cells/1-run_LDA/convert_LDA_to_bins.R"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.no_filt"
outdir=$indir
# indir=$outdir
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt"


for inf in `ls -d $indir/*.Robj`; do
    echo $inf
    bname=$(basename $inf)
    bname=${bname%.*}
    outf=$outdir/$bname.OutObjs.RData
    Rscript $rs $inf $outf&
done
wait
