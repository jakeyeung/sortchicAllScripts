#!/bin/sh
# Jake Yeung
# 0-run.annotate_hiddenDomains_bed.sh
# Run  
# 2019-01-29

# jmark="H3K4me1"
# jmark="H3K4me3"
# jmark="H3K27me3"
jmark="H3K9me3"

inscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/annotate_bed_to_gene_and_distance.R"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output/BM_${jmark}_merged.1000.cutoff/BM_${jmark}_merged.1000.cutoff_analysis.blacklistfilt.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output/BM_${jmark}_merged.1000.cutoff/BM_${jmark}_merged.1000.cutoff_analysis.blacklistfilt.annot.bed"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $wd
Rscript $inscript $inf $outf
