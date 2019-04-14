#!/bin/sh
# Jake Yeung
# run.wig_liftover_bigwig.sh
# Run wig liftover for many files 
# 2019-03-19

rs="/Users/yeung/projects/scchic/scripts/public_data/wig_liftover_bigwig.sh"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
indir="/Users/yeung/data/scchic/public_data/GSE60005_RAW/wigs"
outdir="/Users/yeung/data/scchic/public_data/GSE60005_RAW/bigwigs"

n=0
maxjobs=4

for w in `ls -d $indir/*.wig.gz`; do
		bash $rs $w $outdir&
		if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
			# define maxjobs and n using maxjobsn skeleton
		    wait # wait until all have finished (not optimal, but most times good enough)
		    echo $n wait
		fi
done
wait
echo "Done jobs"
