#!/bin/sh
# Jake Yeung
# 6-make_count_and_RZ_tables.sh
#  
# 2019-12-21

# WRAP UP


# WRAP UP
while [[ `qstat | grep merge | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done


jmem='32G'
jtime='6:00:00'


jsuffix="imputevarfilt.lessstringent.mapq_40"

mapq=40
# binsize=10000  # try other binsizes?
binsize=20000  # try other binsizes?
stepsize=${binsize}
ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.test.py"

# inbed="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_${binsize}.species_drerio.withchr.asinteger.bed"
inbed="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss/zebrafish/gene_tss.CodingOnly.winsize_${binsize}.species_drerio.bed"
inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jsuffix}.remerged_by_marks"
outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.CodingOnly.${jsuffix}.winsize_${binsize}"
[[ ! -d $outdir ]] && mkdir $outdir

for inbami in `ls -d $inmain/*.bam.bai`; do
    inbam=${inbami%.*}
    [[ ! -e $inbam ]] && echo "$inbam not found, continuing" && continue 
    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=$outdir/${bname}.counttables.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    BNAME2=$outdir/${bname}.lhcounts.qsub
    DBASE2=$(dirname "${BNAME2}")
    [[ ! -d $DBASE2 ]] && echo "$DBASE2 not found, exiting" && exit 1

    outf=$outdir/${bname}.countTable.TSS.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu  -N $bname.countTablesTSS
done

