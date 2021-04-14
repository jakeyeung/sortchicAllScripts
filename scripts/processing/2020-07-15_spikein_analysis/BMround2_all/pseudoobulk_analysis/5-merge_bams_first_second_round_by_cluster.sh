#!/bin/sh
# Jake Yeung
# 5-merge_bams_first_second_round_by_cluster.sh
#  
# 2020-11-01

jmem='32G'
jtime='12:00:00'

indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.first_round/all_bams"
indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/split_by_cluster.MAPQ_40.second_round/merged_bams"

# merge common celltypes
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"

jmarks="H3K4me1 H3K4me3 H3K27me3"  # K9me3 separate

jclsts="HSPCs Eryths DCs Bcells Granulocytes NKs"
# jclst="HSPCs"
for jclst in $jclsts; do
    for jmark in $jmarks; do
        BNAME=${outdir}/${jmark}_${jclst}.sbatch_log
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        # merge HSPCs
        infs1=${indir1}/${jmark}*${jclst}*.bam
        infs2=${indir2}/BM_round2_all${jmark}*${jclst}*.bam
        outf=${outdir}/BM_round1_round2_merged_${jmark}_${jclst}.bam
        echo "infs1"
        echo $infs1
        echo "infs2"
        echo $infs2
        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outf $infs1 $infs2; samtools index $outf"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark}_${jclst} --wrap "$cmd"
    done
done

# k9me3 separately
jmark2="H3K9me3"
jclsts2="Eryth Granulocytes HSPCs Lymphoid"
for jclst2 in $jclsts2; do

    BNAME=${outdir}/${jmark2}_${jclst2}.sbatch_log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # merge HSPCs
    infs1=${indir1}/${jmark2}*${jclst2}*.bam
    infs2=${indir2}/BM_round2_all${jmark2}*${jclst2}*.bam
    outf=${outdir}/BM_round1_round2_merged_${jmark2}_${jclst2}.bam
    echo "infs1"
    echo $infs1
    echo "infs2"
    echo $infs2
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools merge $outf $infs1 $infs2; samtools index $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark}_${jclst} --wrap "$cmd"
done
