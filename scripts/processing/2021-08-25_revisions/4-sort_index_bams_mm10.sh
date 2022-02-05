#!/bin/sh
# Jake Yeung
# 4-index_bams_mm10.sh
#  
# 2021-08-25

# WRAP UP
while [[ `squeue -u jyeung | grep BoneMarr | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='2:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/bam_remapped_mm10"
tmpdir="${indir}"

for b in `ls -d $indir/*.bam`; do

    bname=$(basename $b)
    bname=${bname%.*}
    BNAME=${indir}/sbatch_indexbams_log

    bsorted=${indir}/${bname}.sorted.bam

    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools sort -o $bsorted -T $tmpdir $b; samtools index $bsorted"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${bname} --wrap "$cmd"

done


