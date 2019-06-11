#!/bin/sh
# Jake Yeung
# 1-make_bamlist_by_mark.sh
# Create path to bams by mark  
# 2019-05-08

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6-bamlist"
[[ ! -d $outmain ]] && mkdir $outmain

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for mark in $marks; do
    outf=$outmain/JY_${mark}_bamlist_all.out
    [[ -e $outf ]] && echo "$outf already found, exiting" && exit 1
    for indir in $(ls -d $inmain/B6-13W1-BM-${mark}*-merged); do
        echo $indir
        bamdir=$indir/tagged/split_by_bc
        [[ ! -d $bamdir ]] && echo "$bamdir not found, exiting" && exit 1
        for b in $(ls -d $bamdir/*.bam); do
            echo $b >> $outf
        done
    done
done

