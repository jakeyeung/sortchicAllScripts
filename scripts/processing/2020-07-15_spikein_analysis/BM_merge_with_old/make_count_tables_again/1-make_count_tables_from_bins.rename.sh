#!/bin/sh
# Jake Yeung
# 1-make_count_tables_from_peaks.sh
#  
# 2020-11-03

jmem='16G'
jtime='3:00:00'

mapq=40

marks="H3K4me1 H3K9me3 H3K4me3 H3K27me3"

# bamdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/count_tables_from_TSS"
# [[ ! -d $outdir ]] && mkdir $outdir

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_bins/binsize_10000"

for f in `ls -d $inmain/*.txt`; do
    fbase=$(basename $f)
    fbase=`echo $fbase | cut -d"." -f1`
    fnew=${fbase}.binsize_10000.txt
    mv $f $inmain/$fnew
done
