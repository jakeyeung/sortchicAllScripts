#!/bin/sh
# Jake Yeung
# download_dr11_to_dr11_chains.sh
#  
# 2020-04-15

lnk="https://hgdownload.soe.ucsc.edu/goldenPath/danRer10/liftOver/danRer10ToDanRer11.over.chain.gz"
outdir="/hpc/hub_oudenaarden/jyeung/data/databases/chainfiles"

cd $outdir
curl -O $lnk
