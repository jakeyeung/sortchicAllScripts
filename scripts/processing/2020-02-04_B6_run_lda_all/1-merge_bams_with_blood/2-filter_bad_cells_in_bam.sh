#!/bin/sh
# Jake Yeung
# 2-filter_bad_cells_in_bam.sh
#  
# 2020-03-14

# WRAP UP
while [[ `qstat | grep "STDIN" | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='32G'
jtime='6:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
mapq=40

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.withBlood"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.withBlood"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.withBlood.bams_by_cond.MAPQ_${mapq}"
[[ ! -d $outdir ]] && mkdir $outdir

# jmarks="H3K4me3"
jmarks="H3K27me3 H3K9me3"

for jmark in $jmarks; do
    inf=${indir}/"PZ-BM_blood-${jmark}-AllMerged.bam"
    annotinf="${annotdir}/AllMerged_with_blood.${jmark}.txt"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ ! -e $annotinf ]] && echo "$annotinf not found, exiting" && exit 1

    BNAME=${outdir}/split_bam_${jmark}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotinf -outdir $outdir -mapq $mapq --add_chr_prefix" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N split_bam_${jmark}
done

