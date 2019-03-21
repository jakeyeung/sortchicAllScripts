#!/bin/sh
# Jake Yeung
# run.bam_to_bigwig.sh
#  
# 2019-03-20

bs="/Users/yeung/projects/scchic/scripts/scripts_analysis/make_bigwigs/bam_to_bigwig.sh"
indir="/Users/yeung/data/scchic/from_cluster/sorted_bams_2019-03-20"
outdir="/Users/yeung/data/scchic/from_cluster/bigwigs_2019-03-20"

n=0
maxjobs=4
for b in `ls -d $indir/*.bam`; do
	bbase=$(basename $b)
	bbase=${bbase%%.*}
	# echo $bbase
	bout=$outdir/$bbase.bw
	bash $bs $b $bout&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
