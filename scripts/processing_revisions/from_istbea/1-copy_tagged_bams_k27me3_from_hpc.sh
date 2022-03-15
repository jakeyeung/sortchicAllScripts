#!/bin/sh
# Jake Yeung
# 1-copy_tagged_bams_k27me3_from_hpc.sh
#  
# 2022-01-25

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs/tagged_bams/BM_k27me3"
outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/tagged_bams/BM_k27me3"

scp $indir/*SL3* $outdir
