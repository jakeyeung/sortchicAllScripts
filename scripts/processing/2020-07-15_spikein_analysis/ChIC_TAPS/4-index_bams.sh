#!/bin/sh
# Jake Yeung
# 4-index_bams.sh
#  
# 2020-09-25

jmem='16G'
jtime='12:00:00'

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/ChIC_TAPS/tagged_bams/merged_by_mark"
outdir=${indir}

for f in `ls -d $indir/*.bam`; do

    fbase=$(basename $f)
    fbase=${fbase%.*}
    BNAME=${outdir}/${fbase}.log
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools index $f"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${fbase} --wrap "$cmd"
done
