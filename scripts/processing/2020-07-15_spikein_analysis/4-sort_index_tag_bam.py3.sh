#!/bin/sh
# Jake Yeung
# 4-sort_index_tag_bam.sh
#  
# 2019-12-19

# # WRAP UP
# while [[ `qstat | grep PZ |  wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

# WRAP UP
# while [[ `squeue -u jyeung | grep jyeung |  wc -l` > 0 ]]; do
#         echo "sleep for 666 seconds"
#         sleep 666
#     while [[ `squeue -u jyeung | grep jyeung | wc -l` > 0 ]]; do
#         echo "sleep for 66 seconds"
#         sleep 66
#     done
# done

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataK562/raw_demultiplexed"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/raw_demultiplexed"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
# tmpdir=$inmain/tmpdir
tmpdir=$inmain/tmpdir.scmo3
[[ ! -d $tmpdir ]] && mkdir $tmpdir

jmem='64G'
jtime='8:00:00'
ncores=4

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataK562/PZ-K562.tagged_bams"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.py3"
[[ ! -d $outdir ]] && mkdir $outdir

for indir in `ls -d $inmain/PZ*`; do
    bname=$(basename $indir)
    # inbam=$indir/bwaMapped.bam
    inbam=$indir/bwaMapped.py3.bam
    BNAME=$indir/$bname.sort_index_tag.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # sortedbam=$indir/bwaMapped.sorted.bam
    sortedbam=$indir/bwaMapped.py3.sorted.bam
    # outbamtagged=$outdir/${bname}.tagged.bam
    outbamtagged=$outdir/${bname}.py3.tagged.bam
    [[ -e $outbamtagged ]] && echo "$outbamtagged found, continuing" && continue
    # [[ -e $outbam ]] && echo "$outbam found, continuing" && continue
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam; bamtagmultiome.py -method chic --cluster -clusterdir $tmpdir -o $outbamtagged -mem 64 -time 8 $sortedbam" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N sortindextag.$bname
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo3; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam; bamtagmultiome.py -method chic --cluster -clusterdir $tmpdir -o $outbamtagged -mem 64 -time 8 $sortedbam"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done
