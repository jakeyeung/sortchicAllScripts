#!/bin/sh
# Jake Yeung
# 2b-remove_custom_BM_blacklist.sh
# Remove peaks that are highly correlated across the four marks?? 
# 2019-04-21
# generated from "~/projects/scchic/scripts/scripts_analysis/tf_activity_debugging/explore_count_matrix_filter_correlated_peaks.R"

# mark="H3K4me3"
# mark="H3K9me3"
# mark="H3K27me3"
mark="H3K4me1"

suffix="build95_B6_from_traj"

# blacklist="/home/hub_oudenaarden/jyeung/projects/scChiC/data/blacklists/mm10_BM_${mark}_blacklist_QuantileCutoff-0.99.bed.gz"
blacklist="/home/hub_oudenaarden/jyeung/projects/scChiC/data/blacklists/B6_correlated_bins.2019-05-11.bed"   # from 2019-05-12_send_correlated_bins_for_filtering.sh 
inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_${suffix}/merged_across_clusters_${mark}/merged_${mark}.1000.cutoff_analysis.blacklistfilt.bed"
outbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_${suffix}/merged_across_clusters_${mark}/merged_${mark}.1000.cutoff_analysis.blacklistfilt.CorrPeakFilt.bed"

[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
[[ ! -e $blacklist ]] && echo "$blacklist not found, exiting" && exit 1
[[ -e $outbed ]] && echo "$outbed found, exiting" && exit 1

bedtools intersect -v -a $inbed -b $blacklist > $outbed

