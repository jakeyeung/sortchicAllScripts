#!/bin/sh
# Jake Yeung
# 0-make_links_from_archive.sh
#  
# 2019-11-22

indir="/hpc/archive/hub_oudenaarden/seqdata"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
cd $indir

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6LinNeg"
[[ ! -d $outdir ]] && mkdir $outdir

find . -name "*PZ-Bl6-BM-Linneg*" 2>/dev/null | xargs -I {} ln -s {} $outdir/.
