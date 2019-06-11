#!/bin/sh
# Jake Yeung
# 5b-run_make_sitecount_matrix_from_bed.sh
# Avoid SQL nonesense and just load data to memory 
# 2019-02-04

jmem='128G'
jtime='6:00:00'

# jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# jmarks="H3K27me3 H3K9me3"
# jmarks="H3K27me3"
jmarks="H3K4me1"

# threses="0 1000 10000"
threses="0 1000 10000"

n=0
maxjobs=6

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/make_sitecount_matrix_from_bed.R"
# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jscale=0
jcenter=0
byrow=0

prefix="cluster"
suffix="build95_B6_from_traj_CorrPeakFilt"  # suffix for sitecount matrix
suffixCounts="${suffix}.cells_from_bin_analysis"  # MARA analysis now includes exprs matrix with more clels

# fstr="filt.${thres}.promoter"

for thres in $threses; do
    fstr="filt.${thres}"
    for jmark in $jmarks; do
        echo $jmark
        inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene_${prefix}_${suffix}/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.$fstr.long.bed"
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

        outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_${prefix}_${suffixCounts}/${jmark}/mara_input/sitecount_mats"
        [[ ! -d $outdir ]] && mkdir -p $outdir
        outf="${outdir}/hiddenDomains_motevo_merged.closest.long.scale_${jscale}.center_${jcenter}.byrow_${byrow}.bugfix.$fstr.txt"

        BNAME=${outdir}/${fstr}_qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        [[ -e $outf ]] && echo "$outf found, continuing" && continue
        # create the directories beforehand probably
        if [ byrow == 0 ]
        then
            # Rscript $rs $inf $outf --scale $jscale --center $jcenter&
            echo "cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
        else
            # Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow&
            echo "cd $wd; Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
        fi
    done
done
wait
echo "Done script"
