#!/bin/sh
# Jake Yeung
# 2d-merge_counts_across_offsets.sh
# Ran across many offsets. Now merge 
# https://unix.stackexchange.com/questions/60577/concatenate-multiple-files-with-same-header
# 2019-05-06

n=0
maxjobs=40

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

for indir in `ls -d $inmain/B6-13W1-BM-H3K*-merged`; do
    bname=$(basename $indir)
    tagdir=$indir/tagged
    refname=$bname.filtered.bincounts.offset_0.csv
    refinf=$tagdir/$refname
    [[ ! -e $refinf ]] && echo "$refinf not found, exiting" && exit 1
    infs=$(echo $tagdir/$bname.filtered.bincounts.offset_*.csv)
    # infs=$(ls $tagdir/$bname.filtered.counts.offset_*.csv | echo)
    echo $infs
    outf="$tagdir/$bname.filtered.bincounts.slidewin.csv"
    # rm $outf
    [[ -e $outf ]] && echo "$outf  found, exiting" && exit 1
    head -2 $refinf > $outf && tail -q -n +3 $infs >> $outf
    # echo "gzip $outf"
done
wait
