#!/bin/sh
# Jake Yeung
# 1-project_new_samples_on_LDA.H3K27me3.StemCells.sh
#  
# 2019-10-24

jmark="H3K27me3"
jdate="2019-11-17"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_lineage_neg/2-run_LDA/project_new_samples_on_LDA_bin.R"

inlda="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.${jmark}.RData"
inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells_part2/PZ-Bl6-BM-StemCells_matsMerged_H3K27me3_${jdate}.RData"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
[[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-StemCells"
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1
outf="${outdir}/PZ-Bl6-BM-StemCells_matsMerged_${jmark}_${jdate}.RData"

jmem='64G'
jtime='36:00:00'
BNAME=${outdir}/${jmark}_qsub
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

Rscript $rs $inlda $inmat $outf
