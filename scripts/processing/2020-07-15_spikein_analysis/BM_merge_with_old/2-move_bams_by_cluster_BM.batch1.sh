#!/bin/sh
# Jake Yeung
# 2-move_bams_by_cluster_BM.batch1.sh
#  
# 2020-11-04

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/split_by_cluster.MAPQ_40.batch1"
outdir="${inmain}/merged_bam/."

for indir in `ls -d $inmain/H3*`; do
    for inf in `ls -d $indir/*bam*`; do
        ln -s $inf $outdir
    done
done
