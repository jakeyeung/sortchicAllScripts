#!/bin/sh
# Jake Yeung
# 3-merge_peaks_across_clusters.sh
# Merge peaks across clusters 
# 2019-04-15

inbase="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_build95"
[[ ! -d $inbase ]] && echo "$inbase not found, exiting" && exit 1

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
minlength=1000

for mark in $marks; do
    outdir=${inmain}/hd_merged.${mark}.minlength_${minlength}
    [[ ! -d $outdir ]] && mkdir $outdir
    outname="merged.${mark}.minlength_${minlength}.cutoff_analysis.bed"
    outf=${outdir}/${outname}
    [[ -e $outf ]] && echo "$outf found, skipping" && continue
    inmain=${inbase}/hd_clusters.${mark}.minlength_${minlength}
    for indir in $(ls -d ${inmain}/${mark}*); do
        bdir=$(basename $indir)
        inbed="${indir}/${bdir}_analysis.bed"
        [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
        echo "cat $inbed >> $outf"
    done
done
