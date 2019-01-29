#!/bin/sh
# Jake Yeung
# 0-download_fasta_database.sh
# Download fasta database 
# 2019-01-29

outdir="/hpc/hub_oudenaarden/jyeung/data/databases/fasta"
# lnk="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz"

cd $outdir
# curl -O $lnk

# tar -zxvf "chromFa.tar.gz"

# merge fa together
# cat *.fa > "mm10.fa"

# clean up
rm chromFa.tar.gz
rm chr*.fa
