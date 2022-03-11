#!/bin/sh
# Jake Yeung
# 1-run.filter_bed_by_ChIP_peaks.sh
#  
# 2021-06-04

jmem='16G'
jtime='2:00:00'

# compare H3K27me3 vs H3K27me3
jmark="H3K4me2"
# refbed="/hpc/hub_oudenaarden/jyeung/data/public_data/ENCODE/bedfiles/ENCODEpeaks.${jmark}.chr.cleanchr.qval_3.bed"  # bed b
refbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/${jmark}/merged_beds/hiddendomains_output/K562_${jmark}_20181120_allmerged.1000.cutoff/K562_${jmark}_20181120_allmerged.1000.cutoff_analysis.bed"
[[ ! -e $refbed ]] && echo "$refbed not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/${jmark}"

# outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/${jmark}"
[[ ! -d $outmain ]] && echo "outmain $outmain not found, exiting" && exit 1
outdir="${outmain}/Kaya-Okur_vs_hiddendomains_peaks_filtered"
[[ ! -d $outdir ]] && mkdir $outdir

for beda in `ls -d $inmain/*.bed.gz`; do
    bedamain=$(basename $beda)
    bedamain=${bedamain%.*}
    outbed=${outdir}/${bedamain}.ENCODE_filt.bed
    BNAME=${outdir}/${bedamain}.sbatchout
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bedtools intersect -wa -a $beda -b $refbed > $outbed"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=bedtools_${BNAME} --wrap "$cmd"
    # exit 0
done


