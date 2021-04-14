#!/bin/sh
# Jake Yeung
# 1-make_bigwigs_public_data.sh
#  
# 2020-11-22

jmem='16G'
jtime='1:00:00'

indir="/hpc/hub_oudenaarden/Peter/data/K562/published/newbams"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_ENCODE_rerun_from_Buys"
marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.bed"
bs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/bam_to_bigwig_mm10_with_blacklist.sh"
bsize=1000

for mark in $marks; do
    fbase="${mark}_ChIP_dedup_index"
    fname=${fbase}.bam
    inbam=${indir}/${fname}
    outbw=${outdir}/${fbase}.bsize_${bsize}.bw

    BNAME=${outdir}/${fbase}.sbatch
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bash $bs $inbam $outbw ${bsize} $bl"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=bw_${mark} --wrap "$cmd"
done
