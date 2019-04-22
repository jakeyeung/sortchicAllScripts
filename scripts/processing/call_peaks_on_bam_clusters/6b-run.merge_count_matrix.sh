#!/bin/sh
# Jake Yeung
# 6b-run.merge_count_matrix.sh
# Run  
# 2019-04-19

mark="H3K4me1"
tssdist="50000"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/call_peaks_on_bam_clusters/merge_count_matrix.R"

inf1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.1000_build95.withchr.cells_from_bin_analysis/PZ-BM-${mark}.merged.NoCountThres.hiddenDomains.Robj"
inf2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.1000_build95.withchr_${tssdist}.cells_from_bin_analysis/PZ-BM-${mark}.merged.NoCountThres.GeneTSS.Dedup.Robj"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.1000_build95.withchr_${tssdist}.cells_from_bin_analysis/PZ-BM-${mark}.merged.NoCountThres.GeneTSS.Dedup.RbindHiddenDomains.Robj"

[[ ! -e $inf1 ]] && echo "$inf1 not found, exiting" && exit 1
[[ ! -e $inf2 ]] && echo "$inf2 not found, exiting" && exit 1
[[ -e $outf ]] && echo "$outf  found, exiting" && exit 1

Rscript $rs $inf1 $inf2 $outf
