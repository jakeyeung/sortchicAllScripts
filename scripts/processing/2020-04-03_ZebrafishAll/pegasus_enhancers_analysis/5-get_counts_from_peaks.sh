#!/bin/sh
# Jake Yeung
# 4-get_counts_from_peaks.sh
#  
# 2020-04-15

jmem='32G'
jtime='6:00:00'

mapq=40
# dist=10000
dist=20000
ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.test.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inbed="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer11_CNEs_PEGASUS.forliftover.${dist}.nochr.bed"
[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.ConservedEnhancers.winsize_${dist}"
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

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -bedfile $inbed" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.countTable
done
