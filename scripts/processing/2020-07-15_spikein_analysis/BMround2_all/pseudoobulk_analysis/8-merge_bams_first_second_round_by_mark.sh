#!/bin/sh
# Jake Yeung
# 5-merge_bams_first_second_round_by_cluster.sh
#  
# 2020-11-01

jmem='90G'
jtime='32:00:00'

indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.first_round/all_bams"
indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.second_round/merged_bams"

# merge common celltypes
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds.by_mark"
[[ ! -d $outdir ]] && mkdir $outdir

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"  # K9me3 separate
for jmark in $jmarks; do
    BNAME=${outdir}/${jmark}.sbatch_log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # merge HSPCs
    infs1=${indir1}/${jmark}*.bam
    infs2=${indir2}/BM_round2_all${jmark}*.bam
    outf=${outdir}/BM_round1_round2_merged_${jmark}.bam
    echo "infs1"
    echo $infs1
    echo "infs2"
    echo $infs2
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outf $infs1 $infs2; samtools index $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=Merge_${jmark} --wrap "$cmd"
done
