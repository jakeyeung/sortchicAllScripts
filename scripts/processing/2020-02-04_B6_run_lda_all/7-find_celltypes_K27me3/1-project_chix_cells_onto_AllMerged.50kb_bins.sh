#!/bin/sh
# Jake Yeung
# 1-project_new_samples_on_LDA.H3K27me3.StemCells.sh
#  
# 2019-10-24


jmem='24G'
jtime='24:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jmarkref="H3K27me3"
# jmarkref="H3K4me1"
jmarks="K4m1 K27m3 K4m1_K27m3"
# markdbl="K4m1_K27m3"



ldadir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2"
ldaname="lda_outputs.BM_${jmarkref}_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_${jmarkref}_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj"
inlda="${ldadir}/${ldaname}"
[[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_B6_from_chix"
outdir=${indir}/LDA_outputs_projections
[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

inmat="${indir}/H3K27me3_ChIX_padded.rds"
[[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1

BNAME=${outdir}/${jmarkref}_on_AllMerged.qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

outf=${outdir}/${jmarkref}_project_on_BM_AllMerged.RData

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $inmat $outf"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=project_${jmarkref} --wrap "$cmd"

