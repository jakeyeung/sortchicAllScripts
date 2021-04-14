#!/bin/sh
# Jake Yeung
# 3B-reannotate_distance_on_bed_withoutsplitting.sh
#  
# 2020-08-25

jmem='72G'
jtime='4:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/annotate_distance_on_bed.R"
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $wd

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# jmarks="H3K4me3"
# jmarks="H3K27me3"
# jmarks="H3K9me3"
# jmarks="H3K4me1"
jmark="H3K27me3"

n=0
maxjobs=1

    # inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM_dynamic_bins_1000/H3K9me3/motevo_outputs/bed/merged_bed_closestbed_long"
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM_dynamic_bins_1000/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long"
    [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
    fname="motevo_merged.closest.long.bed"
    inf=$inmain/${fname}
        echo $inf
        bname=$(basename $inf)
        # bname=${bname%.*}
        outf=$inmain/${bname}.reannotated.bed
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        [[ -e $outf ]] && echo "$outf found, exiting" && exit 1

        BNAME=${inmain}/${bname}.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=ReannotateDists --wrap "$cmd"
wait
