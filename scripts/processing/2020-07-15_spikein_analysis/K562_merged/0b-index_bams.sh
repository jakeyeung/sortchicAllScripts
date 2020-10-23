#!/bin/sh
# Jake Yeung
# 0b-index_bams.sh
#  
# 2020-09-14

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
cd $indir

jmem='16G'
jtime='3:00:00'

for f in `ls -d $indir/*.bam`; do
    bname=$(basename $f)
    BNAME=${indir}/${bname}.index
    echo $BNAME
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; samtools index $f"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=indx_${bname} --wrap "$cmd"
done
