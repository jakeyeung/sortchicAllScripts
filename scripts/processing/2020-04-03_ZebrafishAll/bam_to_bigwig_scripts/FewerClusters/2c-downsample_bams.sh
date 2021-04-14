#!/bin/sh
# Jake Yeung
# 1b-downsample_bams.sh
#  
# 2020-06-19

# first count reads in all bms

jmem='8G'
jtime='1:00:00'

jprefix="imputevarfilt.lessstringent"
mapq=40
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jprefix}.mapq_${mapq}.FewerClusters.2020-06-17/DedupR1onlyNoAltHits"
outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jprefix}.mapq_${mapq}.FewerClusters.2020-06-17/DedupR1onlyNoAltHits.Downsamp"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do
    outf=${indir}/counts_in_bams.${jmark}.txt
    mincounts=`cut -d"," $outf -f2 | sort -n | head -n 1`
    for inbam in `ls -d $indir/PZ-ChIC-ZF_${jmark}*.bam`; do
        bname=$(basename $inbam)
        bname=${bname%.*}

        BNAME=${outdir}/${bname}.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        outbam=${outdir}/${bname}.bam  # keep same name will make life easier probably
        [[ -e $outbam ]] && echo "$outbam found, continuing" && continue
        frac=$( samtools idxstats $inbam | cut -f3 | awk -v mincounts=$mincounts 'BEGIN {total=0} {total += $1} END {frac=(mincounts*0.99)/total; if (frac > 1) {print 1} else {print frac}}' )
        total=$( samtools idxstats $inbam | cut -f3 | awk -v mincounts=$mincounts 'BEGIN {total=0} {total += $1} END {print total}' )

        echo "Downsampling $bname to $frac calculated by ${mincounts}*0.99 / ${total}"

        echo "frac is $frac,  do subsampling"
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools view -bs $frac $inbam > $outbam; samtools index $outbam"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=ds_${jmark} --wrap "$cmd"
    done
done
