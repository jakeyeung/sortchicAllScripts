#!/bin/sh
# Jake Yeung
# 6-gzip_fastqs.sh
#  
# 2021-06-08

n=0
maxjobs=8

jmem='24G'
jtime='24:00:00'

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/prefetch_outputs/fastqs"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/fastqs"
cd $indir

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/fastqs/gziplogs"
[[ ! -d $outdir ]] && mkdir $outdir

for f in `ls -d $indir/*fastq`; do
    echo $f
    cmd="gzip $f"
    fbase=$(basename $f)
    fbase=${fbase%.*}

    BNAME=${outdir}/${fbase}.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${fbase} --wrap "$cmd"
done
