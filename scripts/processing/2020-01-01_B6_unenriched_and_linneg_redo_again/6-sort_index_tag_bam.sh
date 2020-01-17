#!/bin/sh
# Jake Yeung
# 4-sort_and_index_bam.sh
#  
# 2019-09-04

# # WRAP UP
# while [[ `qstat | wc -l` > 1 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13/raw_demultiplexed"
tmpdir=$inmain/tmpdir
[[ ! -d $tmpdir ]] && mkdir $tmpdir

jmem='64G'
jtime='4:00:00'
ncores=4

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_mergedAll.retag"
# outdir=$inmain/retagged
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams"
[[ ! -d $outdir ]] && mkdir $outdir

for indir in `ls -d $inmain/*BM*`; do
# for indir in `ls -d $inmain/B6-13W1-BM-H3K4me3-1`; do
    bname=$(basename $indir)
    # echo $indir
    # outdir="$indir/tagged"
    inbam=$indir/bwaMapped.bam
    BNAME=$indir/$bname.sort_index_tag.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    sortedbam=$indir/bwaMapped.sorted.bam
    outbamtagged=$outdir/${bname}.retagged.bam
    [[ -e $outbamtagged ]] && echo "$outbamtagged found, continuing" && continue
    # [[ -e $outbam ]] && echo "$outbam found, continuing" && continue
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam; bamtagmultiome.py -method chic --cluster -clusterdir $tmpdir -o $outbamtagged -mem 64 -time 8 $sortedbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N sortindextag.$bname
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam; bamtagmultiome.py -method chic --cluster -o $outbamtagged -mem 32 -time 4 $sortedbam"
done
