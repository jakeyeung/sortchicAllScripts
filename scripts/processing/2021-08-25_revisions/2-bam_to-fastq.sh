#!/bin/sh
# Jake Yeung
# 2-bam_to-fastq.sh
#  
# 2021-08-25

# bam to fastq then remap 

jmem='16G'
jtime='24:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/bam_originals_mm9"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/fastqs_from_bam"

for b in `ls -d $indir/*.bam`; do
    bname=$(basename $b)
    bname=${bname%*}
    echo $bname
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools fastq $b > ${outdir}/${bname}.fastq"

    BNAME=${outdir}/${bname}.sbatch_log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${bname} --wrap "$cmd"
done
