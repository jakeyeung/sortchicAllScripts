#!/bin/sh
# Jake Yeung
# 3-make_merged_bams_from_split_bams.sh
# Splitting bams already filters cells, so make a new merged bam from the split bams (handles bad cells etc) 
# 2020-02-14

jmem='32G'
jtime='6:00:00'

mapq="40"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_${mapq}"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster.MAPQ_${mapq}"
[[ ! -d $outdir ]] && mkdir $outdir

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for mark in $marks; do
    echo $mark
    inbams=${inmain}/${mark}-BM_AllMerged*.sorted.bam
    outbam="${outdir}/${mark}-BM_AllMerged.merged_by_clusters_with_NAs.bam"

    BNAME=$outdir/${mark}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outbam $inbams; samtools index $outbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N Remerged_${mark}
    # echo "samtools merge $outbam $inbams"

done
