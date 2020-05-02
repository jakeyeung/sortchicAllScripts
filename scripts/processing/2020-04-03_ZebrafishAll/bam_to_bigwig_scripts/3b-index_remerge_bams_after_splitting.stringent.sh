#!/bin/sh
# Jake Yeung
# 3-remerge_bams_after_splitting.stringent.sh
#  
# 2020-04-24

jmem='4G'
jtime='1:00:00'

jsuffix="imputevarfilt.lessstringent.mapq_40"

inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jsuffix}.remerged_by_marks"
[[ ! -e $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir=${inmain}

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do

    inbam="${inmain}/${jmark}.imputevarfilt.lessstringent.mapq_40.remerged.bam"
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1

    BNAME=${outdir}/index_remerge_qsub_${jmark}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools index $inbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N remerge_${jmark}

done

