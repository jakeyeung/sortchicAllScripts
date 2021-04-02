#!/bin/sh
# Jake Yeung
# 7-dedup_matrix.sh
#  
# 2020-03-03

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/sitecount_mats_WithIntercept/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.addint_1.marks_merged.WithDupes.WithDupes.txt"
# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/sitecount_mats_WithIntercept/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.addint_1.marks_merged.WithDupes.txt"
# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.marks_merged.txt"
# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.marks_merged.txt"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.marks_merged.WithDupes.txt"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.marks_merged.txt"

awk '!seen[$0]++' ${inf} > ${outf}

