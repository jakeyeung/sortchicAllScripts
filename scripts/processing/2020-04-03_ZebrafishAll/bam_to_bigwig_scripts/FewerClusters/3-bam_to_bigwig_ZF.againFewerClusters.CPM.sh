#!/bin/sh
# Jake Yeung
# 2-bam_to_bigwig.sh
#  
# 2020-01-10

# # WRAP UP
# while [[ `qstat | grep split | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='16G'
jtime='6:00:00'

mapq="40"
jprefix="imputevarfilt.lessstringent"

suffix=".Downsamp"

bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.offset.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_${mapq}.2020-06-17"
# indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jprefix}.mapq_${mapq}.FewerClusters.2020-06-17/DedupR1onlyNoAltHits"
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jprefix}.mapq_${mapq}.FewerClusters.2020-06-17/DedupR1onlyNoAltHits${suffix}"
# outdir="$indir/DedupR1onlyNoAltHits"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/zebrafish/WKM/correlated_bins.chromofixed.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

bsizes="100"
for bsize in $bsizes; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bigwigs_all/ZebrafishWKM.bigwigs_by_fewerclusters.mapq_${mapq}.bsize_${bsize}.2020-06-17.offset.r1only${suffix}"
    [[ ! -d $outdir ]] && mkdir $outdir
    for inbam in `ls -d $indir/*.bam`; do
        echo $inbam
        bname=$(basename $inbam)
        bname=${bname%.*}  # strip extension 
        outbw=${outdir}/${bname}.${bsize}.cleaned.CPM.bw
        [[ -e $outbw ]] && echo "$outbw found, continuing" && continue

        BNAME=$outdir/$bname.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw $bsize $bl"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
        # echo $outdir
        # echo $outbw
    done
done
