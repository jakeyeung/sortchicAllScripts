#!/bin/sh
# Jake Yeung
# 4b-run_make_exprs_matrix_from_LDA_peaks.sh
#  
# 2020-08-25

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/make_merged_exprs_mat_for_MARA.args.R"

jmark="H3K4me1"

indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZFWKM_peaks/lda_outputs.${jmark}.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.binarize.FALSE"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/${jmark}/mara_input"
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

inf="${indir}/ldaOut.${jmark}.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.Robj"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

# inf="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZFWKM_peaks/lda_outputs.H3K4me1.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.binarize.FALSE/ldaOut.H3K4me1.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.Robj"

bname=$(basename $inf)
bname=${bname%.*}

outdir="${outmain}/count_mats_peaks_norm_merged"
[[ ! -d $outdir ]] && mkdir $outdir
outf="${outdir}"/${bname}.txt

BNAME=${outdir}/${bname}.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
