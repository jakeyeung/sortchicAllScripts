#!/bin/sh
# Jake Yeung
# 4b-run_make_exprs_matrix_from_LDA_peaks.sh
#  
# 2020-08-25

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/make_merged_exprs_mat_for_MARA.args.R"

keepNbins=0  # keeps all bins
# jmark="H3K4me1"
# jmark="H3K4me3"
# jmark="H3K27me3"
jmarks="H3K4me1 H3K4me3 H3K27me3"

for jmark in $jmarks; do

    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_topbins"
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged3_Peaks"
    [[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

    jsuffix="lda_outputs.count_mat_from_hiddendomains.${jmark}.filtNAcells_topbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.${jmark}.filtNAcells_topbins.K-30.Robj"
    inf=${inmain}/${jsuffix}
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    bname=$(basename $inf)
    bname=${bname%.*}

    outdir="${outmain}/count_mats_peaks_norm"
    [[ ! -d $outdir ]] && mkdir $outdir
    outf="${outdir}"/${bname}.keepNbins_${keepNbins}.txt
    [[ -e $outf ]] && echo "$outf found, exiting" && exit 1

    BNAME=${outdir}/${bname}.keepNbins_${keepNbins}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf -keepNbins $keepNbins --AddChr" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1

done
