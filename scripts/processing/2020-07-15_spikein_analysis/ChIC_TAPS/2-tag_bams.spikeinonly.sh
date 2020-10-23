#!/bin/sh
# Jake Yeung
# 4-sort_index_tag_bam.sh
#  
# 2019-12-19

# # WRAP UP
# while [[ `squeue -u jyeung | grep Bwa |  wc -l` > 0 ]]; do
#         echo "sleep for 66 seconds"
#         sleep 66
# done


jmem='64G'
jtime='8:00:00'
ncores=4

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams/spikeins_only"
[[ ! -d $outdir ]] && mkdir $outdir

tmpdir=$outdir/tmpdir
[[ ! -d $tmpdir ]] && mkdir $tmpdir

indir="/hpc/hub_oudenaarden/avo/chictaps/VAN4973/lambda"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS"
# [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

# for indir in `ls -d $indir/*.lambda.sorted.bam`; do
    # bname=$(basename $indir)
    # inbam=$indir/bwaMapped.bam
    for inbam in `ls -d $indir/*.lambda.sorted.bam`; do
        bname=$(basename $inbam)
        bname=${bname%%.*}  # clip all dots out
        BNAME=${outdir}/$bname.sort_index_tag.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        sortedbam=$outdir/bwaMapped.${bname}.sorted.bam
        outbamtagged=$outdir/${bname}.contigfixed.tagged.bam
        [[ -e $outbamtagged ]] && echo "$outbamtagged found, continuing" && continue
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam; bamtagmultiome.py -method chic --cluster -clusterdir $tmpdir -o $outbamtagged -mem 64 -time 8 $sortedbam"
        # cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; bamtagmultiome.py -method chic --cluster -clusterdir $tmpdir -o $outbamtagged -mem 64 -time 8 $sortedbam"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
    done
# done
