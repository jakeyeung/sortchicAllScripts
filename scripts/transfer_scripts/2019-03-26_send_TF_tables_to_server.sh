#!/bin/sh
# Jake Yeung
# 2019-03-26_send_TF_tables_to_server.sh
#  
# 2019-03-26

# indir="/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.H3K4me1_activities_only.RData"
indir="/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.H3K4me1_activities_only.with_count_mats.RData"
outdir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook"

scp $indir $outdir

