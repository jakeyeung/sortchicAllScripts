#!/bin/sh
# Jake Yeung
# 3-remerge_after_split.sh
#  
# 2020-03-14

# WRAP UP
while [[ `qstat | grep split_bam | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='32G'
jtime='6:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.withBlood.bams_by_cond.MAPQ_40"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.withBlood.bams_remerged.MAPQ_40"
[[ ! -d $outdir ]] && mkdir $outdir

# jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# jmarks="H3K4me1 H3K4me3"
jmarks="H3K27me3 H3K9me3"

for jmark in ${jmarks}; do

    outbam=${outdir}/PZ-BM_blood-${jmark}-AllMerged.remerged.bam
    inbams=${indir}/PZ*${jmark}*.bam
    echo $inbams

    # check inbams exist
    for inbam in $inbams; do
        [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    done

    BNAME=${outdir}/remerge_${jmark}.bam
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outbam $inbams; samtools index $outbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N remerge_${jmark} -m beas -M j.yeung@hubrecht.eu

done
