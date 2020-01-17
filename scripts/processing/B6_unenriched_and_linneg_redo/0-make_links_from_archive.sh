#!/bin/sh
# Jake Yeung
# 0-make_links_from_archive.sh
#  
# 2019-11-22

indir="/hpc/archive/hub_oudenaarden/seqdata"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
cd $indir

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6"
[[ ! -d $outdir ]] && mkdir $outdir

find . -name "B6-13W1-BM-*" 2>/dev/null | xargs -I {} ln -s {} $outdir/.
