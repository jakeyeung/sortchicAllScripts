#!/bin/sh
# Jake Yeung
# 1-merge_marks_with_blood.sh
#  
# 2020-03-14

jmem='64G'
jtime='12:00:00'

bmdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final"
bldir="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/raw_data/3-analysis_intestines.tagged_bams.all_merged"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.withBlood"
[[ ! -d $outdir ]] && mkdir $outdir
jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do
    b1="${bmdir}/${jmark}-BM_AllMerged.bam"
    b2=${bldir}/"PZ_blood_${jmark}_2020-02-29.bam"

    outbase="PZ-BM_blood-${jmark}-AllMerged"
    outbam=${outdir}/${outbase}.bam

    BNAME=${outdir}/${outbase}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools merge $outbam $b1 $b2; samtools index $outbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
done
