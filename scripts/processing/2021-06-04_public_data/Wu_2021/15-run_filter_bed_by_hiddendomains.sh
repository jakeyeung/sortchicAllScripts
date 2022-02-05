#!/bin/sh
# Jake Yeung
# 1-run.filter_bed_by_ChIP_peaks.sh
#  
# 2021-06-04

jmem='16G'
jtime='2:00:00'

refbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams_r1_notrim_bowtie2/hiddendomains_output/SRR12638101_1.sorted.2500.cutoff/SRR12638101_1.sorted.2500.cutoff_analysis.bed"
[[ ! -e $refbed ]] && echo "$refbed not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/${jmark}"
# inmain=""
beda1="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/GSM4780540_K27me3_H1_r1.fragments.HG38.nochr.tsv"
# beda2="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/bams_r1_notrim_bowtie2/hiddendomains_output/SRR12638101_1.sorted.2500.cutoff/SRR12638101_1.sorted.2500.cutoff_analysis.bed"
beda2="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data/prefetch_outputs/SRR12638101/counts_output_r1_notrim_bowtie2/SRR12638101_1.sorted_dupcounts.bed"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/peak_filtered_beds_again"
[[ ! -d $outdir ]] && mkdir $outdir

for beda in ${beda1} ${beda2}; do
    bedamain=$(basename $beda)
    bedamain=${bedamain%.*}
    outbed=${outdir}/${bedamain}.hiddendomains_filt.bed
    BNAME=${outdir}/${bedamain}.sbatchout
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bedtools intersect -wa -a $beda -b $refbed > $outbed"
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=bedtools_${BNAME} --wrap "$cmd"
    # exit 0
done
