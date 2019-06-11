#!/bin/sh
# Jake Yeung
# 8-make_TSS_peak_filtered_list.sh
# Take input from LDA to get fasta later 
# 2019-04-19

# mark="H3K4me3"
# mark="H3K27me3"
# mark="H3K9me3"
mark="H3K4me1"
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="$wd/scripts/processing/call_peaks_on_bam_clusters/make_bed_from_count_mat.R"
suffix="CorrPeakFilt"
prefix="build95_B6_from_traj"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.1000_${prefix}.cells_from_bin_analysis/${suffix}"
fname="B6-${mark}.merged.NoCountThres.hiddenDomains.${suffix}.Robj"
inf=$indir/$fname

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_${prefix}/merged_across_clusters_${mark}"
fnameout="merged_${mark}.1000.cutoff_analysis.blacklistfilt.${suffix}.fromCountMat.bed"
outf=$outdir/$fnameout

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

cd $wd; Rscript $rs $inf $outf

