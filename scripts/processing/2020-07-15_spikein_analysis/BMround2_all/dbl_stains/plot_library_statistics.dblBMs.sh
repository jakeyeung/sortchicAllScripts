#!/bin/sh
# Jake Yeung
# plot_library_statistics.sh
#  
# 2020-08-17

jmem='16G'
jtime='3:00:00'

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/merged_bams"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_dbl_stains/dbl_stains"

outmain=${indir}/libstats_by_plate
[[ ! -d $outmain ]] && mkdir $outmain

for inbam in `ls -d $indir/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}


    outdir=${outmain}/${bname}
    [[ ! -d $outdir ]] && mkdir $outdir

    BNAME=${outdir}/${bname}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; libraryStatistics.py $inbam -o $outdir"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${bname} --wrap "$cmd"
done

