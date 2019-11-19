#!/bin/sh
# Jake Yeung
# 6-make_RZ_counts.sh
#  
# 2019-09-04

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_HD_0/raw_demultiplexed"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_HD_0/RZcounts"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='8G'
jtime='2:00:00'

# for indir in `ls -d $inmain/PZ*`; do
for indir in `ls -d $inmain/PZ-ChIC-ZF*`; do
    bname=$(basename $indir)
    inbam=$indir/tagged/bwaMapped.sorted.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, skipping" && continue
    outf=$outdir/${bname}.RZ_counts.csv

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py $inbam -sampleTags SM -featureTags RZ -o $outf --dedup --filterXA -minMQ 40" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname.RZcounts
done

