#!/bin/sh
# Jake Yeung
# 2d-merge_counts_across_offsets.sh
# Ran across many offsets. Now merge 
# https://unix.stackexchange.com/questions/60577/concatenate-multiple-files-with-same-header
# 2019-05-06

n=0
maxjobs=40

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

for indir in `ls -d $inmain/B6-13W1-BM-H3K*-merged`; do
    bname=$(basename $indir)
    tagdir=$indir/tagged
    outf="$tagdir/$bname.filtered.bincounts.slidewin.csv"
    [[ ! -e $outf ]] && echo "$outf not found, exiting" && exit 1
    gzip $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
