#!/bin/sh
# Jake Yeung
# 1-copy_bams_to_istbea.sh
#  
# 2022-07-18

inf="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds"

outf="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/bams/notH3K27me3"

scp $inf/*.bam $outf
