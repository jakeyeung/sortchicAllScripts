#!/bin/sh
# Jake Yeung
# 6-fastqdump.sh
#  
# 2021-06-09

jmem='16G'
jtime='1:00:00'

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/fastq_outputs"
[[ ! -d $outmain ]] && mkdir $outmain

for d in `ls -d $inmain/SRR*`; do
    dbase=$(basename $d)
    outdir=${outmain}/${dbase}
    [[ ! -d $outdir ]] && mkdir $outdir
    
    inf=${d}/${dbase}.sra

    BNAME=${outdir}/${dname}.sbatch.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd="fastq-dump --outdir ${outdir} --split-files $inf --gzip"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${dbase} --wrap "$cmd"
done
