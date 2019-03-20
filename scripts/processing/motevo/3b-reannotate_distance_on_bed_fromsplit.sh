#!/bin/sh
# Jake Yeung
# 3b-reannotate_distance_on_bed_fromsplit.sh
# Annotate distances from split 
# 2019-03-07
#  

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/annotate_distance_on_bed.R"
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $wd

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# jmarks="H3K4me1"
jmarks="H3K27me3"

n=0
maxjobs=12

for jmark in $jmarks; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$jmark/motevo_outputs/bed/merged_bed_closestbed_long/split"
    for inf in `ls -d $inmain/motevo*`; do
        echo $inf
        bname=$(basename $inf)
        outf=$inmain/$bname.bed
        # outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$jmark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        [[ -e $outf ]] && echo "$outf found, exiting" && exit 1
        Rscript $rs $inf $outf&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        	# define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done
wait
