#!/bin/sh
# Jake Yeung
# 3-merge_peaks_across_clusters.sh
# Merge peaks across clusters 
# 2019-04-15

suffix="build95_B6_from_traj"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_${suffix}"

marks="H3K4me1"
distfilt=1000

for mark in $marks; do
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_${suffix}/merged_across_clusters_${mark}"
    [[ ! -d $outdir ]] && mkdir $outdir
    outname="merged_${mark}.${distfilt}.cutoff_analysis.blacklistfilt.bed"
    outf=${outdir}/${outname}
    [[ -e $outf ]] && echo "$outf found, skipping" && continue
    for indir in $(ls -d ${inmain}/${mark}_cluster_*); do
        bdir=$(basename $indir)
        inbed="${indir}/${bdir}_analysis.blacklistfilt.bed"
        [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
        cat $inbed >> $outf
    done
done
