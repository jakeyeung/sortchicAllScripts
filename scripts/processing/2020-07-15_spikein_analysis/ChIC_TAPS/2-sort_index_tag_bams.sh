#!/bin/sh
# Jake Yeung
# 4-sort_index_tag_bam.sh
#  
# 2019-12-19

# WRAP UP
while [[ `squeue -u jyeung | grep Bwa |  wc -l` > 0 ]]; do
        echo "sleep for 66 seconds"
        sleep 66
done

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
tmpdir=$inmain/tmpdir
[[ ! -d $tmpdir ]] && mkdir $tmpdir

jmem='64G'
jtime='8:00:00'
ncores=4

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams"
[[ ! -d $outdir ]] && mkdir $outdir

for indir in `ls -d $inmain/CG*`; do
    bname=$(basename $indir)
    inbam=$indir/bwaMapped.bam
    BNAME=$indir/$bname.sort_index_tag.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    sortedbam=$indir/bwaMapped.sorted.bam
    outbamtagged=$outdir/${bname}.contigfixed.tagged.bam
    [[ -e $outbamtagged ]] && echo "$outbamtagged found, continuing" && continue
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam; bamtagmultiome.py -method chic --cluster -clusterdir $tmpdir -o $outbamtagged -mem 64 -time 8 $sortedbam"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done
