#!/bin/sh
# Jake Yeung
# 3-plot_umaps_sweep_params.sh
# Plot umap output of LDA 
# 2019-05-08

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/2-run_LDA_pipeline_on_binned_matrix_build95_B6/plot_umap.R"

thres="0.995"
ncores=16
nnmin=30
nnmax=62
nnby=1
mindist=0.4
k="50"

# jbin="TRUE"; kvec="25_30_35_50"
# jbin="TRUE"; kvec="25_30_40_50"
# jbin="FALSE"; kvec="30_40_50"
# jbin="TRUE"; kvec="15_17_20_23"
# jbin='FALSE'
jbin='FALSE'

# marks="H3K9me3"
# marks="H3K9me3"
marks="H3K27me3 H3K4me1 H3K4me3 H3K9me3"
# marks="H3K27me3 H3K9me3"

for mark in $marks; do
    echo $mark
    # indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BinFiltCellFilt/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.${jbin}.no_filt"
    indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.${jbin}.no_filt"
    for indat in `ls -d $indir/lda_out_meanfilt.B6_${mark}_pcutoff*.Robj`; do
        # indat="$indir/lda_out_meanfilt.B6_${mark}_pcutoff_0.CountThres0.K-${kvec}.Robj"
        echo $indat


        bname=$(basename $indat)
        bname=${bname%.*}
        outpdf="${indir}/${bname}_mindist_${mindist}_umap_plots_winbugfix.pdf"
        outrdata="${indir}/${bname}_mindist_${mindist}_mindist_processed_lda.Rdata"

        [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
        # [[ -e $outpdf ]] && echo "$outpdf found, continuing" && exit 1
        # [[ -e $outrdata ]] && echo "$outrdata found, exiting" && exit 1

        # echo $outpdf
        # echo "cd $wd; Rscript $rs $indat $outpdf $outrdata $mark $thres $ncores $nnmin $nnmax $nnby $mindist $k"
        cd $wd; Rscript $rs $indat $outpdf $outrdata $mark $thres $ncores $nnmin $nnmax $nnby $mindist $k
    done
done

