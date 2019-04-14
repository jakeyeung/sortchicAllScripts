#!/bin/sh
# Jake Yeung
# run.variance_over_pseudotime_plots_parallel_for_server.sh
#  
# 2019-04-04

jmem='42G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/scripts_analysis/build95_primetime_from_server/variance_over_pseudotime_plots_parallel_for_server.OneMark.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/plot_outputs/spatial_analysis_again"

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
for mark in $marks; do
    outpdf=${outmain}/variance_over_pseudotime_plots.${mark}.pdf
    BNAME=$outmain/nohup_${mark}.pdf
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo "cd $wd; Rscript $rs $mark $outpdf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err 
done
