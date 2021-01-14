#!/bin/sh
# Jake Yeung
# 3b-reannotate_distance_on_bed_fromsplit.sh
# Annotate distances from split 
# 2019-03-07
#  

jmem='8G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/annotate_distance_on_bed.R"
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $wd

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# jmarks="H3K4me3"
# jmarks="H3K27me3"
# jmarks="H3K9me3"
jmarks="H3K4me1 H3K4me3 H3K27me3"

n=0
maxjobs=1

for jmark in $jmarks; do
    # inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged_Peaks_1000/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/split"
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged3_Peaks/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/split"
    for inf in `ls -d $inmain/*_split.*`; do
        echo $inf
        bname=$(basename $inf)
        # bname=${bname%.*}
        outf=$inmain/${bname}.reannotated.bed
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        [[ -e $outf ]] && echo "$outf found, exiting" && exit 1

        BNAME=${inmain}/${bname}.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N Annot_${jmark}
        # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf
    done
done
wait
