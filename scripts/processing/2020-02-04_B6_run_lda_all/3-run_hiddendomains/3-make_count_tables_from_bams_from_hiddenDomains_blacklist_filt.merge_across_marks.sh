#!/bin/sh
# Jake Yeung
# 2-make_count_tables_from_bams_from_hiddenDomains_blacklist_filt.sh
#  
# 2020-02-14

# WRAP UP
while [[ `qstat | grep "Remerged_H" | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='96G'
jtime='12:00:00'
mapq=40

minlength=1000

ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.test.py"  # allow blacklist as option 
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

marks="H3K4me1 H3K27me3 H3K4me3"

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
[[ ! -e ${bl} ]] && echo "${bl} not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/count_tables_from_hiddenDomains_marks_merged"
[[ ! -d $outdir ]] && mkdir $outdir

# bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final"
bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster"

hdout="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/hd_merged.marks_merged/marks_merged.ntopic_30.binskeep_250.txt"
[[ ! -e $hdout ]] && echo "$hdout not found, exiting" && exit 1

for mark in ${marks}; do
    echo $mark
    inbam=${bamdir}/"${mark}-BM_AllMerged.merged_by_clusters_no_NAs.bam"
    [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1

    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outftab=${outdir}/${bname}.txt
    [[ -e $outftab ]] && echo "$outftab not found, continuing" && continue

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps --filterXA -minMQ $mapq $inbam -o $outftab -sampleTags SM -joinedFeatureTags reference_name  -binTag DS --dedup -bed $hdout -blacklist $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N MakeTable_${bname}
        
    # exit 0
done
