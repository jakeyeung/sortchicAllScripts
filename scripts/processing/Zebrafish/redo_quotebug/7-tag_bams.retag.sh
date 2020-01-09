#!/bin/sh
# Jake Yeung
# 4-tag_bams.sh
#  
# 2019-09-04

# WRAP UP

# sleep 7200  # 2 hour wait

# ggGfastarefunzipped="Danio_rerio.GRCz11.dna.primary_assembly.unzipped.fa"
# fastaref="/hpc/hub_oudenaarden/group_references/ensembl/98/danio_rerio/Danio_rerio.GRCz11.dna.primary_assembly.unzipped.fa"
fastaref="/hpc/hub_oudenaarden/group_references/ensembl/98/danio_rerio/Danio_rerio.GRCz11.dna.primary_assembly.bgzipformat.fa.gz"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all/raw_demultiplexed.first_try"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

tmpdir="$inmain/tmpdir"
[[ ! -d $tmpdir ]] && mkdir $tmpdir

jmem='8G'
jtime='12:00:00'

for indir in `ls -d $inmain/PZ-ChIC-ZF*`; do
    bname=$(basename $indir)
    echo $indir
    tmpprefix=$tmpdir/$bname
    outdir="$indir/tagged.retagged"
    [[ ! -d $outdir ]] && mkdir $outdir
    # skip if outbam is not empty
    inbam=$indir/tagged/${bname}.bwaMapped.tagged.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, continuing" && continue
    outbam=$outdir/${bname}.bwaMapped.tagged.retagged.bam  # only works if using bamtagmultiome
    [[ -e $outbam ]] && echo "$outbam found, continuing" && continue
    BNAME=$indir/$bname.tag.retag.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamtagmultiome.py -method chic --cluster -o $outbam -mem 8 $inbam -ref $fastaref"  | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N tag.$bname

    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; cd $tmpdir; bamtagmultiome.py -method chic -o $outbam $inbam -ref $fastaref&

    # exit 0  # test one to see if identical
    # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamtagmultiome.py -method chic -o $outbam -mem 8 $inbam -ref $fastaref
    # exit 0
     # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; universalBamTagger -tmpprefix $tmpprefix --chic --ftag -moleculeRadius 5 -o $outdir $inbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname
done
wait
