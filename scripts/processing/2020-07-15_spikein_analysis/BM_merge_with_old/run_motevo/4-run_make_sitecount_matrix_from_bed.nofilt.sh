#!/bin/sh
# Jake Yeung
# 5b-run_make_sitecount_matrix_from_bed.sh
# Avoid SQL nonesense and just load data to memory 
# 2019-02-04

jmem='72G'
jtime='6:00:00'

jmarks="H3K4me1 H3K4me3 H3K27me3"

n=0
maxjobs=6

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/motevo_scripts/lib/make_sitecount_matrix_from_bed.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# dont scale or center, but scale and center the expression matrix?
jscale=0
jcenter=0
byrow=0

jsuffix="BM-AllMerged3_Peaks"

for jmark in $jmarks; do
    echo $jmark
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_${jsuffix}/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_${jsuffix}/${jmark}/mara_input/sitecount_mats"
    [[ ! -d $outdir ]] && mkdir -p $outdir

    outf="${outdir}/hiddenDomains_motevo_merged.closest.long.scale_${jscale}.center_${jcenter}.byrow_${byrow}.${jmark}.txt"

    BNAME=${outdir}/${jmark}_qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # create the directories beforehand probably
    echo $byrow
    
    if [ $byrow -eq 0 ]
    then
        # Rscript $rs $inf $outf --scale $jscale --center $jcenter&
        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N SCMat_${jmark}
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter"
    else
        # Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow&
        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -N SCMatnorm_${jmark}
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow"
    fi
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=SCMat_${jmark} --wrap "$cmd"
done
wait
echo "Done script"
