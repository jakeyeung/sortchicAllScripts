#!/bin/sh
# Jake Yeung
# download_dr7_to_dr10_chains.sh
#  
# 2020-04-15

lnk="https://hgdownload.soe.ucsc.edu/goldenPath/danRer7/liftOver/danRer7ToDanRer10.over.chain.gz"
outdir="/hpc/hub_oudenaarden/jyeung/data/databases/chainfiles"

cd $outdir
curl -O $lnk
