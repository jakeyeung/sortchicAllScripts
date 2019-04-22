#!/bin/sh
# Jake Yeung
# 5-run.select_best_gene.sh
#  
# 2019-04-19

mark="H3K4me1"
tssdist=60000
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/call_peaks_on_bam_clusters/select_best_gene.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.1000_build95.withchr_${tssdist}.cells_from_bin_analysis/PZ-BM-${mark}.merged.NoCountThres.GeneTSS.Robj"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.1000_build95.withchr_${tssdist}.cells_from_bin_analysis/PZ-BM-${mark}.merged.NoCountThres.GeneTSS.Dedup.Robj"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ -e $outf ]] && echo "$outf found, exiting" && exit 1
Rscript $rs $inf $outf


