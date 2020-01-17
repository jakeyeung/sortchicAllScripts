#!/bin/sh
# Jake Yeung
# 4-tag_bams.sh
#  
# 2019-09-04

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_mergedAll.retag"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_mergedAll.retag.universalBamTagger"
[[ ! -d $outdir ]] && mkdir $outdir

tmpdir="$outdir/tmpdir"
[[ ! -d $tmpdir ]] && mkdir $tmpdir

jmem='32G'
jtime='3:00:00'

for inbam in `ls -d $inmain/*.retagged.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    echo $bname
    tmpprefix=$tmpdir/$bname
    # skip if outbam is not empty
    outbam=$outdir/${bname}.retagged.univeraslBamTagger.bam
    [[ -e $outbam ]] && echo "$outbam found, continuing" && continue
    # [[ -s $outbam ]] && echo "$outbam not empty, continuing" && continue
    BNAME=$indir/$bname.tag.debug.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamtagmultiome.py -method chic --cluster -o $outbam -mem 8 -time 3 $inbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N tag.$bname
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; cd $outdir; bamtagmultiome.py -method chic --cluster -o $outbam -mem 32 -time 3 $inbam
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; universalBamTagger.py -tmpprefix $tmpprefix --chic --ftag -moleculeRadius 5 -o $outdir $inbam&
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamtagmultiome.py -method chic --cluster -o $outbam -mem 8 -time 3 $inbam"
done
wait

