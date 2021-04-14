#!/bin/sh
# Jake Yeung
# 2-run_LDA_merged_with_chix.sh
#  
# 2020-09-03
jmem='8G'
jtime='48:00:00'


rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

topics=30
bname="H3K27me3_AllMerged_and_ChIX_merged"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_B6_from_chix/H3K27me3_AllMerged_and_ChIX_merged.rds"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_B6_from_chix/LDA_AllMerged_with_ChIX"
[[ ! -d $outdir ]] && mkdir $outdir

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname"

sbatch --time=$jtime --mem-per-cpu=$jmem --output=${outdir}/${bname}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
