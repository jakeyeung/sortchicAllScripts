#!/bin/sh
# Jake Yeung
# 2-select_best_gene.sh
#  
# 2019-06-24

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/call_peaks_on_bam_clusters/select_best_gene.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3700/countTables_genebodies"
outdir=$indir
# prefixs="mm_K36me3-K27me3 mm_H3K36me3-H3K9me3 mm_H3K4me1 mm_H3K27me3 mm_H3K9me3"
# for prefix in $prefixs; do
for inf in `ls -d $indir/*2000.RData`; do
    bname=$(basename $inf)
    bname=${bname%.*}
    outf=$outdir/$bname.dedup.RData
    Rscript $rs $inf $outf&
done
wait
