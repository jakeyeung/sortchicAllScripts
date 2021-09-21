#!/bin/sh
# Jake Yeung
# 0-download_refseq_TES.sh
#  
# 2020-08-10

link="https://ccg.epfl.ch/mga/mm10/refseq/MmRefseqTes.sga.gz"

outdir="/hpc/hub_oudenaarden/jyeung/data/databases/refseq"

cd $outdir

curl -O $link
