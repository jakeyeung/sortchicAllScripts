#!/bin/sh
# Jake Yeung
# 9-sort_index_bams.sh
#  
# 2021-06-10

jmem='16G'
jtime='4:00:00'

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/bams"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/bams_demux_bugfix"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams"
tmpdir="/hpc/hub_oudenaarden/jyeung/data/tmp"
ncores=4

# for d in `ls -d $indir/SRR*`; do 
for inbam in `ls -d $indir/*.bam`; do 
    # inbam=${indir}/${dbase}/${dbase}.bam
    # [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
    dbase=$(basename $inbam)
    dbase=${dbase%.*}
    sortedbam=${indir}/${dbase}.sorted.bam

    # BNAME=${indir}/${dbase}/${dbase}.sbatch.log
    BNAME=${indir}/${dbase}.sort_index.sbatch.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; samtools sort -T $tmpdir -@ $ncores $inbam > $sortedbam; samtools index $sortedbam"
    echo $cmd
    # exit 0
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${dbase} --wrap "$cmd"

done
