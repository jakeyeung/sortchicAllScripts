#!/bin/sh
# Jake Yeung
# 4-tag_bams.sh
#  
# 2019-09-04

# ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_unenriched_and_linneg_redo/universalBamTagger.addtempdir.sh"
# [[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6LinNeg/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

tmpdir="$inmain/tmpdir"
[[ ! -d $tmpdir ]] && mkdir $tmpdir

jmem='8G'
jtime='12:00:00'

for indir in `ls -d $inmain/*BM*`; do
    bname=$(basename $indir)
    echo $indir
    tmpprefix=$tmpdir/$bname
    outdir="$indir/tagged"
    # skip if outbam is not empty
    # overwrite so dont check empty
    # outbam=$outdir/bwaMapped.sorted.bam
    # [[ -s $outbam ]] && echo "$outbam not empty, continuing" && continue

    inbam=$indir/bwaMapped.sorted.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, continuing" && continue
    BNAME=$indir/$bname.tag.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; universalBamTagger -tmpprefix $tmpprefix --chic --ftag -moleculeRadius 5 -o $outdir $inbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamtagmultiome.py -method chic --cluster -o $outdir -mem 8 $inbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname
done

