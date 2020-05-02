#!/bin/sh
# Jake Yeung
# 3-remerge_bams_after_splitting.stringent.sh
#  
# 2020-04-24

jmem='64G'
jtime='6:00:00'

jsuffix="imputevarfilt.lessstringent.mapq_40"
inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jsuffix}"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jsuffix}.remerged_by_marks"
# [[ ! -d $outdir ]] && mkdir $outdir

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

cd $inmain

for jmark in $jmarks; do
    inbams=${inmain}/PZ-ChIC-ZF*${jmark}*sorted.bam
    outbam="${outdir}/${jmark}.${jsuffix}.remerged.bam"

    echo $jmark
    echo $inbams
    echo $outbam

    BNAME=${outdir}/remerge_qsub_${jmark}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outbam $inbams" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N remerge_${jmark}
done
