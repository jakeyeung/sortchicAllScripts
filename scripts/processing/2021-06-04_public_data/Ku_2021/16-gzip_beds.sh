#!/bin/sh
# Jake Yeung
# 16-gzip_beds.sh
#  
# 2021-06-15

jmem='2G'
jtime='1:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/counts_output/blacklist_filt"

for f in `ls -d $indir/*.bed`; do
    bname=$(basename $f)
    bname=${bname%.*}
    BNAME=${indir}/${bname}.gzip.sbatch.log

    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd="gzip $f"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done
