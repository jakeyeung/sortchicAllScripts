#!/bin/sh
# Jake Yeung
# 9-run.summarize_histone_marks_compare_metacell_lib.sh
# Summarize a folder
# 2019-01-13

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/summarize_histone_marks_compare_metacell_lib.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all"
mcdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_dropbox/FIG4_BM"
jdist=1000
jmean=1
jmin=100
jmax=500000
jtops="15_20_25_30_35"
jthres=0.96  # for GREAT
binarize="TRUE"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_downstream/ldaAnalysisHiddenDomains_${jdist}/lda_outputs.meanfilt_${jmean}.cellmin_${jmin}.cellmax_${jmax}.binarize.${binarize}"
[[ -d $outdir ]] && echo "$outdir found, exiting for safety" && exit 1
[[ ! -d $outdir ]] && mkdir -p $outdir
echo $outdir
cd $wd; Rscript $rs $indir $outdir $mcdir $jdist $jmean $jmin $jmax $jtops $jthres $binarize
