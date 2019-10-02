#!/bin/sh
# Jake Yeung
# 3-downstream_LDA_sweep_params.sh
#  
# 2019-06-18

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/2-run_LDA_pipeline_on_binned_matrix_build95_B6/plot_umap_do_not_add_chr.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

thres="0.995"
ncores=16
nnmin=30
nnmax=62
nnby=1
mindist=0.4
k="50"

jbin='TRUE'

prefix="PZ-Bl6-BM-Linneg"
marks="H3K4me3"

for mark in $marks; do
    echo $mark
    indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_${prefix}/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.${jbin}.no_filt"
    # for indat in `ls -d $indir/lda_out_meanfilt.B6_${mark}_pcutoff*.Robj`; do
    for indat in `ls -d $indir/lda_out_meanfilt.PZ-Bl6-BM-Linneg_${mark}_binfilt_cellfilt.CountThres0.*.Robj`; do
        echo $indat
        bname=$(basename $indat)
        bname=${bname%.*}
        outpdf="${indir}/${bname}_mindist_${mindist}_umap_plots_winbugfix.pdf"
        outrdata="${indir}/${bname}_mindist_${mindist}_mindist_processed_lda.Rdata"
        [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
        cd $wd; Rscript $rs $indat $outpdf $outrdata $mark $thres $ncores $nnmin $nnmax $nnby $mindist $k
    done
done
