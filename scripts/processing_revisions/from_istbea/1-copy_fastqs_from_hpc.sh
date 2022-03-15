#!/bin/sh
# Jake Yeung
# 1-copy_fastqs_from_hpc.sh
#  
# 2022-01-21

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/raw_data"
outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC"

# scp -r $indir $outdir
rsync -avrL --partial --partial-dir=$tmpdir --copy-links $indir $outdir
