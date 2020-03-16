#!/bin/sh
# Jake Yeung
# 5b-run_make_sitecount_matrix_from_bed.sh
# Avoid SQL nonesense and just load data to memory 
# 2019-02-04

jmem='128G'
jtime='6:00:00'

jmarks="H3K27me3"

n=0
maxjobs=6

wd="/home/jyeung/projects/scchic"
rs="/home/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/make_sitecount_matrix_from_bed.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# dont scale or center, but scale and center the expression matrix?
jscale=0
jcenter=0
byrow=0

jsuffix="BM-AllMerged_Peaks_1000"

for jmark in $jmarks; do
    echo $jmark
    inf="/home/jyeung/hpc/scChiC/tfbs_output_cluster_${jsuffix}/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    outdir="/home/jyeung/hpc/scChiC/mara_analysis_${jsuffix}/${jmark}/mara_input/sitecount_mats"
    [[ ! -d $outdir ]] && mkdir -p $outdir

    outf="${outdir}/hiddenDomains_motevo_merged.closest.long.scale_${jscale}.center_${jcenter}.byrow_${byrow}.${jmark}.txt"

    BNAME=${outdir}/${fstr}_qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # create the directories beforehand probably
    echo $byrow
    
    if [ $byrow -eq 0 ]
    then
        # Rscript $rs $inf $outf --scale $jscale --center $jcenter&
        # echo "cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N SCMat_${jmark}
        cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter
    else
        # Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow&
        # echo "cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N SCMatnorm_${jmark}
        cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow
    fi
done
wait
echo "Done script"
