#!/bin/sh
# Jake Yeung
# 1-copy_bams_to_istbea.sh
#  
# 2022-07-18

inf="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster/bams_merged_by_cluster"

outf="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/bams"

scp $inf/*.bam $outf
