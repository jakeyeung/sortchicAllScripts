#!/bin/sh
# Jake Yeung
# filter_bc_matches_B6.sh
# Filter out only cells with BC matches 
# 2019-05-06

bcfile="/home/hub_oudenaarden/jyeung/projects/scChiC/data/barcode_summaries/barcodes/maya_384NLA.bc"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_B6/parsed2"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_B6/parsed2_filt"

[[ ! -d $outdir ]] && mkdir $outdir

n=0
maxjobs=32

for f in `ls -d $indir/*.out.gz`; do
    fbase=$(basename $f)
    outf=$outdir/$fbase
    zcat $f | fgrep -f $bcfile | gzip > $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
