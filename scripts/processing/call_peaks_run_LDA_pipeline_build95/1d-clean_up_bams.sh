#!/bin/sh
# Jake Yeung
# 1c-clean_up_bams.sh
# Clean up bams after adding prefix
# https://wiki.library.ucsf.edu/display/StateLab/Changing+the+header+of+a+BAM+file 
# now we double bam, remove duplicate
# 2019-03-25

inf1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95/BM_H3K4me1_merged.bam"

inf2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95/BM_H3K27me3_merged.bam"

inf3="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95/BM_H3K9me3_merged.bam"

inf4="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95/BM_H3K4me3_merged.bam"

rm $inf1 $inf2 $inf3 $inf4
