#!/bin/sh
# Jake Yeung
# 1-copy_all_zebrafish_fastqs.sh
#  
# 2020-04-03

inmain="/hpc/archive/hub_oudenaarden/seqdata"
cd $inmain
outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/raw_data"
[[ ! -d $outdir ]] && mkdir $outdir

find . -name "*PZ-ChIC-ZF*.fastq.gz" 2>/dev/null | xargs -I {} cp {} $outdir/.

