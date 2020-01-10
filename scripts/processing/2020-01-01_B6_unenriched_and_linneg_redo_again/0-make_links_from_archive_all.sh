#!/bin/sh
# Jake Yeung
# 0-make_links_from_archive_all.sh
#  
# 2019-12-13

# copy everything

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13"
[[ ! -d $outdir ]] && mkdir $outdir

indir="/hpc/archive/hub_oudenaarden/seqdata"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
cd $indir

find . -name "B6-13W1-BM-*.fastq.gz" 2>/dev/null | xargs -I {} cp {} $outdir/.  # the original wildtype data
find . -name "*PZ-Bl6-BM-Linneg*.fastq.gz" 2>/dev/null | xargs -I {} cp {} $outdir/.  # linneg for all except K4me1 (K4me1 was bad?)
find . -name "*ChIC-B6BMSC*.fastq.gz"  2>/dev/null | xargs -I {} cp {} $outdir/.  # K4me3, K27me3, K9me3 stem cells 
find . -name "PZ-ChIC-Bl6-BM-*.fastq.gz" 2>/dev/null | xargs -I {} cp {} $outdir/.  # H3K4me1 new antibody 
