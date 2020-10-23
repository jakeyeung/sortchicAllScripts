#!/bin/sh
# Jake Yeung
# 0-filter_CTCF_peaks_by_chromo.sh
#  
# 2020-08-17

chromos="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX"

inbed="/hpc/hub_oudenaarden/jyeung/data/databases/beds/CTCF_motifs/mm10_bonemarrow-all_CTCF-chip_optimal-IDRpeaks_ENCFF806PDR.bed"
outbed="/hpc/hub_oudenaarden/jyeung/data/databases/beds/CTCF_motifs/mm10_bonemarrow-all_CTCF-chip_optimal-IDRpeaks_ENCFF806PDR.chromofilt.bed"

[[ -e $outbed ]] && echo "$outbed not found, exiting" && exit 1

for c in $chromos; do
    cstr="^$c"
    grep $cstr $inbed >> $outbed
done

# chromos=`for c in `seq 19`; do echo $c; done`
# echo $chromos

