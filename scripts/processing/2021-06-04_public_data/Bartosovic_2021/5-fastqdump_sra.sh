#!/bin/sh
# Jake Yeung
# 5-fastqdump_sra.sh
#  
# 2021-06-07

jmem='16G'
jtime='1:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/prefetch_outputs"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/prefetch_outputs/fastqs"

for d in `ls -d $indir/SRR*`; do
    echo $d
    dname=$(basename $d)
    inf=${d}/${dname}.sra
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    BNAME=${outdir}/${dname}.sbatch.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd="fastq-dump --outdir ${outdir} --split-files $inf"
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${dname} --wrap "$cmd"
    # exit 0
done
