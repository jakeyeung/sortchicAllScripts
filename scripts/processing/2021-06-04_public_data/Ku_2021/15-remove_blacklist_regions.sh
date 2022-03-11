#!/bin/sh
# Jake Yeung
# 15-remove_blacklist_regions.sh
#  
# 2021-06-15

jmem='16G'
jtime='1:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output/blacklist_filt"
blacklist="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.nochr.bed"
[[ ! -e $blacklist ]] && echo "$blacklist not found, exiting" && exit 1

for inbed in `ls -d $indir/*.bed.gz`; do
    fbase=$(basename $inbed)
    fbase=${fbase%%.*}
    outbase="${fbase}.sorted_dupcounts.blacklist_filt.bed"

    BNAME=${outdir}/${outbase}.sbatch.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    outbed=${outdir}/${outbase}
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bedtools subtract -a $inbed -b $blacklist -A > $outbed"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${fbase} --wrap "$cmd"
done
