#!/bin/sh
# Jake Yeung
# parse_cellbarcode_and_dinucleotide.sh
#  
# 2019-04-16

n=0
maxjobs=8
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_VAN2979"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_VAN2979/parsed2"
outdir=$indir/parsed2

[[ ! -d $outdir ]] && mkdir $outdir

# PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14_L002_R1_001
# H3K4me1.BM.m1.S9.AH3VGVBGX9.CAGAATAT
for f in `ls -d $indir/*.out`; do
    bname=$(basename $f)
    outname=${bname%%.*}.parsed.out
    outf=$outdir/$outname
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    awk -F $'\t' -v bname=$bname 'BEGIN {OFS = FS} {split(bname,a,"-"); split(bname,b,"_");  print a[4]"."a[2]"."a[3]"."b[3]"."b[2]"."substr($1,0,8), substr($1,0,8), substr($1,9,9)}' $f > $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
