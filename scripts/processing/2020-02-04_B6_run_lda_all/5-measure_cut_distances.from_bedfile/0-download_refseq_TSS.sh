#!/bin/sh
# Jake Yeung
# 0-download_refseq_TSS.sh
#  
# 2020-08-06

# Download refseq TSS from Bucher

# link="https://ccg.epfl.ch/mga/mm10/refseq/refseq.html"
link="https://ccg.epfl.ch/mga/mm10/refseq/MmRefseqTss.sga.gz"
outdir="/hpc/hub_oudenaarden/jyeung/data/databases/refseq"

cd $outdir

curl -O $link
