#!/bin/sh
# Jake Yeung
# 4-run.plotHeatmap.sh
#  
# 2020-02-24

# WRAP UP
while [[ `qstat | grep compMatrix | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='8G'
jtime='1:00:00'

jmark="H3K27me3"
bsize="bsize_100"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/bigwig_outputs/merged_bams.deeptools_outputs.tss.dist_5000.ctype.sorted.${jmark}.${bsize}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir=${indir}

for inf in `ls -d $indir/*.tab.gz`; do
    bname=$(basename $inf)
    bname=${bname%.*}
    jstr=$(echo $bname | cut -d"." -f2,3)
    BNAME=${outdir}/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outf=${outdir}/plotHeatmap.${bname}.png
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    echo $jstr
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; plotHeatmap -m $inf -out $outf --colorMap YlGnBu --heatmapHeight 15 --sortRegions descend --sortUsingSamples 1 --averageTypeSummaryPlot median" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N ${jstr}
done

