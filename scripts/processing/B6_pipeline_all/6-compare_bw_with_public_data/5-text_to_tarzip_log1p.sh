#!/bin/sh
# Jake Yeung
# 5-text_to_tarzip.sh
# Zip up files for transfer 
# 2019-03-29

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk"
# indir="merged_softlinks_textfile"
indir="merged_softlinks_textfile_log1p_build95_B6"
outf="$inmain/merged_softlinks_textfile_log1p_build95_B6.tar.gz"

[[ ! -d $inmain/$indir ]] && echo "$inmain/$indir not found, exiting" && exit 1

cd $inmain
tar -zcvf $outf $indir
