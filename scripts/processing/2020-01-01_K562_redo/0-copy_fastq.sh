#!/bin/sh
# Jake Yeung
# 0-copy_fastq.sh
#  
# 2019-12-24

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataK562"
[[ ! -d $outdir ]] && mkdir $outdir

inmain="/hpc/archive/hub_oudenaarden/seqdata"
cd $inmain

find . -iname "*PZ-K562-G1-H3*-rep*.fastq.gz" 2>/dev/null | xargs -I {} cp --no-clobber {} $outdir/.  # the original wildtype data
find . -iname "*PZ-ChIC-K562*-no*.fastq.gz" 2>/dev/null | xargs -I {} cp --no-clobber {} $outdir/.  # neg controls

