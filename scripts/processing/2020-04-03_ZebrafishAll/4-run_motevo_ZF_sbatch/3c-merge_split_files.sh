#!/bin/sh
# Jake Yeung
# 3d-merge_split_files.sh
# After calculating distances, remerge the split files 
# 2019-03-07

# marks="H3K4me1 H3K4me3 H3K27me3"
# marks="H3K27me3 H3K9me3"
# marks="H3K4me1"
# marks="H3K4me3"
# marks="H3K27me3"
# marks="H3K9me3"
# marks="H3K27me3"
marks="H3K4me3"

for mark in $marks; do
    # inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene_${suffix}/${mark}/motevo_outputs/bed/merged_bed_closestbed_long"
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_ZF-AllMerged2_Peaks_1000/${mark}/motevo_outputs/bed/merged_bed_closestbed_long"
    # inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged_Peaks_1000/${mark}/motevo_outputs/bed/merged_bed_closestbed_long"
    [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
    indir=$inmain/split/reannotated
    [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
    outf=$inmain/motevo_merged.closest.dist.long.bed
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    cat $indir/*.bed > $outf
done

