#!/bin/sh
# Jake Yeung
# 3-normalize_bigwig_to_input.sh
# Normalize chipseq bigwig to input 
# 2019-03-21

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
binsize=100

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bcells/bigwigs/renamed"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bcells/bigwigs/renamed_relative_to_input"

# Make everything relative to input

# marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
marks="H3K4me1 H3K4me3 H3K27me3"  # no H3K9me3 for Bcells
ctypes="MatBcell ProB HSC"

n=0
maxjobs=8

for ctype in $ctypes; do
    inputbw=$indir/"input_${ctype}_rep1.bw"  # always rep1 by programmatically set
    [[ ! -e $inputbw ]] && echo "$inputbw not found, exiting" && exit 1
    for mark in $marks; do
        if [[ $ctype == "HSC" && $mark == "H3K27me3" ]]
        then
        	continue
        fi
        for inbw in `ls -d $indir/${mark}_${ctype}_*.bw`; do
            inbase=$(basename $inbw)
            inbase=${inbase%.*}
            rep=$(echo $inbase | cut -d"_" -f3)
            outf="$outdir/${mark}_${ctype}_${rep}_reltoinput.bw"
            bigwigCompare -b1 $inbw -b2 $inputbw -bs=$binsize -o $outf&
            if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
            	# define maxjobs and n using maxjobsn skeleton
                wait # wait until all have finished (not optimal, but most times good enough)
                echo $n wait
            fi
        done
    done
done

# inbw="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bcells/bigwigs/GSM1463445_Mat_Bcell_H3K27me3_rep1.bw"
# inputbw="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bcells/bigwigs/GSM1463447_Mat_Bcell_input.bw"



