#!/bin/sh
# Jake Yeung
# run.liftover_bigwig.sh
# Run wig liftover for many files 
# 2019-03-20

rs="/Users/yeung/projects/scchic/scripts/public_data/bedgraph_liftover_bigwig.sh"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
indir="/Users/yeung/data/scchic/public_data/Neutrophils/Gong_GenesAndDev_2017/bedGraph_mm9"
outdir="/Users/yeung/data/scchic/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10"

n=0
maxjobs=4

for w in `ls -d $indir/*.bedGraph.gz`; do
		echo "bash $rs $w $outdir&"
		bash $rs $w $outdir&
		if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
			# define maxjobs and n using maxjobsn skeleton
		    wait # wait until all have finished (not optimal, but most times good enough)
		    echo $n wait
		fi
done
wait
echo "Done jobs"
