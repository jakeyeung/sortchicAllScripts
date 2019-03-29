#!/bin/sh
# Jake Yeung
# 4d-bedgraph_to_bigwig.sh
# bedGraphToBigWig 
# 2019-03-21

chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.txt"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bedgraphs_2019-03-21/log1p"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_2019-03-20/log1p_bigwigs"

[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

n=0
maxjobs=15
for inbed in `ls -d $inmain/*.bedgraph`; do
    # echo $inbed
    bbase=$(basename $inbed)
    bbase=${bbase%.*}
    outbed=$outmain/$bbase.log1p.bw
    bedGraphToBigWig $inbed $chromsizes $outbed&
if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
	# define maxjobs and n using maxjobsn skeleton
    wait # wait until all have finished (not optimal, but most times good enough)
    echo $n wait
fi
done
wait
