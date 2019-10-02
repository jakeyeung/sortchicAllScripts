#!/bin/sh
# Jake Yeung
# 4-tag_bams.sh
#  
# 2019-09-04

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3697/raw_demultiplexed"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/fastq/raw_demultiplexed"

jmem='8G'
jtime='12:00:00'

for indir in `ls -d $inmain/PZ*`; do
    bname=$(basename $indir)
    echo $indir
    outdir="$indir/tagged"
    inbam=$indir/bwaMapped.sorted.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, continuing" && continue
    BNAME=$indir/$bname.tag.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate singlecellmultiomicsenv2; universalBamTagger.py --chic --ftag -moleculeRadius 5 -o $outdir $inbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname
    # echo "universalBamTagger.py --chic --ftag -moleculeRadius 5 -o $outdir $inbam"
done

