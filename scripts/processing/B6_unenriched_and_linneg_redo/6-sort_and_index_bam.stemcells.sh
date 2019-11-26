#!/bin/sh
# Jake Yeung
# 4-sort_and_index_bam.sh
#  
# 2019-09-04

while [[ `qstat | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6StemCells/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
tmpdir=$inmain/tmpdir
[[ ! -d $tmpdir ]] && mkdir $tmpdir

jmem='8G'
jtime='1:00:00'
ncores=4

for indir in `ls -d $inmain/*BM*`; do
    bname=$(basename $indir)
    echo $indir
    outdir="$indir/tagged"
    inbam=$indir/bwaMapped.bam
    BNAME=$indir/$bname.sort.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    outbam=$indir/bwaMapped.sorted.bam
    # [[ -e $outbam ]] && echo "$outbam found, continuing" && continue
    echo "samtools sort -T $tmpdir -@ $ncores $inbam > $outbam; samtools index $outbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N $bname.sortindex
done
