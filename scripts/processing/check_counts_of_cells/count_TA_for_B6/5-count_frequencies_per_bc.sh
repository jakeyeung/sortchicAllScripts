#!/bin/sh
# Jake Yeung
# 5-count_frequencies_per_bc.sh
# Count frequencies per bc 
# 2019-05-06

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_B6/parsed2_filt"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_B6/parsed2_filt_counts"

[[ ! -d $outdir ]] && mkdir $outdir

n=0
maxjobs=32

for f in `ls -d $indir/*.parsed.out.gz`; do 
    # echo $f
    fname=$(basename $f)
    outf=$outdir/$fname
    zcat $f | sort | uniq -c | gzip > $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
