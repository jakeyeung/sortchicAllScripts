#!/bin/sh
# Jake Yeung
# 3-normalize_bigwig_to_input.sh
# Normalize chipseq bigwig to input 
# 2019-03-21

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
binsize=100

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10/renamed"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10/renamed_relative_to_input"
# Make everything relative to input

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
ctypes="Neu Pro"

n=0
maxjobs=4
for ctype in $ctypes; do
    inputbw=$indir/"input_${ctype}.bw"
    [[ ! -e $inputbw ]] && echo "$inputbw not found, exiting" && exit 1
    for mark in $marks; do
        inbw=$indir/${mark}_${ctype}.bw
        outf="$outdir/${mark}_${ctype}_reltoinput.bw"
        echo "bigwigCompare -b1 $inbw -b2 $inputbw -bs=$binsize -o $outf&"
        bigwigCompare -b1 $inbw -b2 $inputbw -bs=$binsize -o $outf&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
            # define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done

