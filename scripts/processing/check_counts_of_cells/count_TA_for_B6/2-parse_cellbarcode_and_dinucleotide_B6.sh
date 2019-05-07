#!/bin/sh
# Jake Yeung
# parse_cellbarcode_and_dinucleotide.sh
#  
# 2019-04-16

n=0
maxjobs=16
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_B6"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir=$indir/parsed2

[[ ! -d $outdir ]] && mkdir $outdir

for f in `ls -d $indir/B6*R1*.out`; do
    bname=$(basename $f)
    outname=${bname%%.*}.parsed.out
    outf=$outdir/$outname
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    awk -F $'\t' -v bname=$bname 'BEGIN {OFS = FS} {split(bname,a,"-"); split(a[5],b,"_");  print a[4]"."a[1]".m1."b[1]"."b[2]"."substr($1,0,8), substr($1,0,8), substr($1,9,9)}' $f > $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
