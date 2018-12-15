#!/bin/sh
# Jake Yeung
# 2-bin_bam_by_windows.sh
# Bin bam by windows rather than macs2 
# 2018-12-14

inbam="/hpc/hub_oudenaarden/jyeung/data/histone-mods/VAN2979_BM_m1_H3K4me3/merged_bams/VAN2979_BM_m1_H3K4me3-merged.sorted.bam"
genomewin="/hpc/hub_oudenaarden/jyeung/data/databases/genomebins/mm10.winsize.5000.stepsize.1000.bed"
# gfile="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.txt"

bname=$(basename $inbam)
bname=${bname%%.*}
outf="/hpc/hub_oudenaarden/jyeung/data/histone-mods/VAN2979_BM_m1_H3K4me3/coverage_out/$bname.bed"

bedtools coverage -a $genomewin -b $inbam > $outf
