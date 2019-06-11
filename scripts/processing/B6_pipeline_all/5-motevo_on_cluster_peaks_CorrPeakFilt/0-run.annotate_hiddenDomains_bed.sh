#!/bin/sh
# Jake Yeung
# 0-run.annotate_hiddenDomains_bed.sh
# Run  
# 2019-01-29

# jmark="H3K4me1"
# jmark="H3K4me3"
# jmark="H3K27me3"
# jmark="H3K9me3"

n=0
# maxjobs=8
maxjobs=1
# marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# marks="H3K4me1"
# marks="H3K4me3"
marks="H3K27me3 H3K9me3"
# marks="H3K9me3"

suffix="build95_B6"
inscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/annotate_bed_to_gene_and_distance.R"
[[ ! -e $inscript ]] && echo "$inscript not found, exiting" && exit 1
for jmark in $marks; do

    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_${suffix}/merged_across_clusters_${jmark}/merged_${jmark}.1000.cutoff_analysis.blacklistfilt.CorrPeakFilt.fromCountMat.bed"
    outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_${suffix}/merged_across_clusters_${jmark}/merged_${jmark}.1000.cutoff_analysis.blacklistfilt.CorrPeakFilt.fromCountMat.annot.bed"

    echo $inf

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ -e $outf ]] && echo "$outf found, exiting" && exit 1

    wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
    cd $wd
    Rscript $inscript $inf $outf&
done
wait
