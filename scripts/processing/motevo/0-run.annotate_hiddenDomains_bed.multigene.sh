#!/bin/sh
# Jake Yeung
# 0-run.annotate_hiddenDomains_bed.multigene.sh
# Run  
# 2019-01-29

inscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/annotate_bed_to_gene_and_distance.R"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output/BM_H3K4me1_merged.1000.cutoff/BM_H3K4me1_merged.1000.cutoff_analysis.blacklistfilt.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output/BM_H3K4me1_merged.1000.cutoff/BM_H3K4me1_merged.1000.cutoff_analysis.blacklistfilt.annot.multigene.debug.bed"

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $wd
Rscript $inscript $inf $outf --dist 50000 --multi
