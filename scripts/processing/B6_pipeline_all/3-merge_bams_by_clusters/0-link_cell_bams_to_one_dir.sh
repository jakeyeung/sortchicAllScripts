#!/bin/sh
# Jake Yeung
# 0-link_cell_bams_to_one_dir.sh
# Link cell bams to one directory 
# 2019-05-09

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"
outdir="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6-split_by_bc"
[[ ! -d $outdir ]] && mkdir $outdir

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for d in `ls -d $inmain/B6*`; do
    bamdir="$d/tagged/split_by_bc"
    [[ ! -d $bamdir ]] && echo "$bamdir not found, exiting" && exit 1
    for b in `ls -d $bamdir/*bam*`; do 
        ln -s $b $outdir/. 
    done
done
