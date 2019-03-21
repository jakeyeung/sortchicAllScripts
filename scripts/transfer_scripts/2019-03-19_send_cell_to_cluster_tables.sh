#!/bin/sh
# Jake Yeung
# 2019-03-19_send_cell_to_cluster_tables.sh
# Send cell-to-cluster tables to server so we can merge bams 
# 2019-03-19

indir="/Users/yeung/data/scchic/tables/bamlist_for_merging"
outdir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook"

scp -r $indir $outdir
