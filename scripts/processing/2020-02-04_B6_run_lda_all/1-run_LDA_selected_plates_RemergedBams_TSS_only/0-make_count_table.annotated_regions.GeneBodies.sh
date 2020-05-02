#!/bin/sh
# Jake Yeung
# 1-bam_to_count_table_TSS.sh
#  
# 2019-06-23

jmem='32G'
jtime='24:00:00'

mapq=40
# dist=1
prefixs="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

# outsuffix=".AnnotatedGeneRegionsWithRepeats"
outsuffix=".AnnotatedGeneRegionsWithGeneBodies"

ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.WithBlacklist.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster.MAPQ_40"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

# inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/ENCODE/all_winsize_${dist}.bed"
# inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/ENCODE/promsenhs_winsize_${dist}.2020-03-29.bed"
inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/ENCODE/genes.2020-03-30.bed"
[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1

bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS${outsuffix}"
[[ ! -d $outdir ]] && mkdir $outdir

for prefix in $prefixs; do
    inbam="${inmain}/${prefix}*.bam"
    outf=${outdir}/${prefix}.countTableTSS.mapq_${mapq}.gene_body.blfiltered.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$outdir/${prefix}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed -blacklist $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu  -N $prefix.countTablesTSS

done
wait
