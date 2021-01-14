#!/bin/sh
# Jake Yeung
# 8-add_chr_to_hidden_domains_outputs.sh
#  
# 2020-11-05

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"
outdir=${indir}

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for mark in $marks; do
    fname="hd_merged.${mark}.FromR.maxcount_40_60_R/merged.${mark}.cutoff_analysis.merged.nochr.bed"
    outname="hd_merged.${mark}.FromR.maxcount_40_60_R/merged.${mark}.cutoff_analysis.merged.withchr2.bed"
    inf=${indir}/${fname}
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    outf=${outdir}/${outname}
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    sed 's/^/chr/g' $inf > $outf
done
