#!/bin/sh
# Jake Yeung
# 0-rename_bigwigs_IdoAmit.sh
# Rename bigwigs 
# 2019-03-26
# Unzip only H3K4me1 and H3K4me3 data

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/bigwig_mm10"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/renamed"

[[ ! -d $outdir ]] && mkdir $outdir

# name format: Mark_Celltype.bw

marks="H3K4me1 H3K4me3"

for mark in $marks; do
    for b in `ls -d $indir/*${mark}*.bw`; do
        bname=$(basename $b)
        bname=${bname%%.*}  # remove .ucsc.bigWig suffix
        markcell=$(echo $bname | cut -d"_" -f2)
        # cell=$(echo $markcell | cut -d"${mark}" -f2)
        cell=$(echo $markcell | sed -e s/"${mark}"//g)
        # echo $mark
        # echo $cell
        # echo $markcell
        # echo $bname
        newname=${mark}_${cell}.bw
        ln -s $b $outdir/$newname
    done
done

