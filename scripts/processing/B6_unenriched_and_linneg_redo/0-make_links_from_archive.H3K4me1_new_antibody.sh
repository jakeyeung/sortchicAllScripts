#!/bin/sh
# Jake Yeung
# 0-make_links_from_archive.sh
#  
# 2019-11-22

indir="/hpc/archive/hub_oudenaarden/seqdata"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6"
[[ ! -d $outdir ]] && mkdir $outdir

cd $indir
# find . -name "PZ-ChIC-Bl6-BM-H3K4me1-*.fastq.gz" 2>/dev/null | xargs -I {} ln -s {} $outdir/.
# find . -name "PZ-ChIC-Bl6-BM-H3K4me1-*.fastq.gz" 2>/dev/null | xargs -I {} ln -s {} $outdir/.
find . -name "PZ-ChIC-Bl6-BM-H3K4me1-*.fastq.gz" 2>/dev/null | xargs -I {} cp {} $outdir/.
