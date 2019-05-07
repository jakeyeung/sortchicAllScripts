#!/bin/sh
# Jake Yeung
# zip_files.sh
# Zip files because slow 
# 2019-05-06

n=0
maxjobs=16
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_B6/parsed2"

for f in `ls -d $indir/*.out`; do
    gzip $f&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
