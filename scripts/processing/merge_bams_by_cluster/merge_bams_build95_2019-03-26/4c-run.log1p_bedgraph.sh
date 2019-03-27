#!/bin/sh
# Jake Yeung
# 4c-run.log1p_bedgraph.sh
# Log1p transform
# 2019-03-21

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bedgraphs_2019-03-21"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bedgraphs_2019-03-21/log1p"

n=0
maxjobs=8
for inbed in `ls -d $inmain/*.bedgraph`; do
    bbase=$(basename $inbed)
    outbed=$outmain/$bbase
    awk -v OFS='\t' '{print $1,$2,$3,log($4+1)/log(2)}' $inbed | LC_COLLATE=C sort -k1,1 -k2,2n > $outbed&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done


