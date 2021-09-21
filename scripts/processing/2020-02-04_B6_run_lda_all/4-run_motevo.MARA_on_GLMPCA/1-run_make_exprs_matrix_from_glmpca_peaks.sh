#!/bin/sh
# Jake Yeung
# 4b-run_make_exprs_matrix_from_LDA_peaks.sh
#  
# 2020-08-25

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/make_merged_exprs_mat_from_GLMPCA_for_MARA.args.R"
# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/make_merged_exprs_mat_for_MARA.args.R"

# jmark="H3K27me3"
# jmark="H3K4me1"
jmark="H3K4me3"

# indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZFWKM_peaks/lda_outputs.${jmark}.imputevarfilt.lessstringent.mapq_40.countTable.HiddenDomains.NewCountFilters.K-30.binarize.FALSE"
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.Faster.poi.ZF_AllMerged.imputevarfilt.lessstringent.LDAfromPeaks"


# jmark <- "H3K4me1"

# indir="/hpc/hub_oudenaarden/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.Faster.poisKeepBestPlates2.LDAfromPeaks"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/GLMPCA_outputs.Faster.poisKeepBestPlates2.LDAfromPeaks"

fname="PZ_${jmark}.AllMerged.KeepBestPlates2.LDAfromPeaks.GLMPCA_var_correction.mergebinsize_100.binskeep_250.covar_cell.var.within.sum.norm.log2.CenteredAndScaled.penalty_1.5.winsorize_FALSE.2020-08-26.RData"

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/${jmark}/mara_input"
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

inf=${indir}/${fname}
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

bname="BM_${jmark}.BM_AllMerged.VarCorrection.Poisson.GLMPCA"

outdir="${outmain}/count_mats_peaks_from_GLMPCA"
[[ ! -d $outdir ]] && mkdir $outdir
outf="${outdir}"/${bname}.txt
outpdf="${outdir}"/${bname}.pdf

BNAME=${outdir}/${bname}.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

# echo ${BNAME}
# echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -outpdf $outpdf" 
echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf -outfile $outf -outpdf $outpdf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu
