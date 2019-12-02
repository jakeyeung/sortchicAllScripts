#!/bin/sh
# Jake Yeung
# 7-make_count_tables.sh
#  
# 2019-09-04


inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_HD_0/raw_demultiplexed"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_HD_0/countTables"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='9G'
jtime='2:00:00'

mapq=40
stepsize=10000
binsize=50000


# for indir in `ls -d $inmain/PZ*`; do
for indir in `ls -d $inmain/PZ-ChIC-ZF*`; do  # zebrafish only 
    bname=$(basename $indir)
    inbam=$indir/tagged/bwaMapped.sorted.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, skipping" && continue
    outf=$outdir/${bname}.countTable.csv

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.countTable
done

