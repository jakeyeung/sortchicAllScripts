#!/bin/sh
# Jake Yeung
# 2019-03-20_copy_merged_bams.sh
# Copy merged bams to visualize 
# 2019-03-20

inf1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_2019-03-20/H3K4me1_cluster_4.bam"
inf2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_2019-03-20/H3K4me3_cluster_5.bam"
inf3="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_2019-03-20/H3K27me3_cluster_5.bam"
inf4="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_2019-03-20/H3K9me3_cluster_2.bam"

infs="$inf1 $inf2 $inf3 $inf4"

for inf in $infs; do
		scp $inf $outdir
done
