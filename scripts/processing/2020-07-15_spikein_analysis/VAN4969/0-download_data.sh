#!/bin/sh
# Jake Yeung
# 0-download_data.sh
#  
# 2020-08-10

# VAN4969.tar.gpg
# link="https://ncie01.op.umcutrecht.nl/index.php/s/Qa5rp46R9WtGTEy"
link="https://ncie01.op.umcutrecht.nl/s/Qa5rp46R9WtGTEy/download"
outdir="/hpc/hub_oudenaarden/seqdata/VAN4969"

cd $outdir

curl -O $link
