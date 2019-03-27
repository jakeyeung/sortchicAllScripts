#!/bin/sh
# Jake Yeung
# 1c-add_chromo_prefix.sh
# Change header on bams
# 2019-03-25
# 

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95/add_chromo_prefix"

for bam_name in $(ls -d $indir/*_merged.bam); do
    echo $bam_name
    bname=$(basename ${bam_name})
    outf=$outdir/$bname
    samtools view -H $bam_name | sed "s/SN\:/SN\:chr/" | samtools reheader - ${bam_name} > $outf
done


