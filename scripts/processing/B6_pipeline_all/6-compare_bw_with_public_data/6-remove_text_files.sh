#!/bin/sh
# Jake Yeung
# 6-remove_text_files.sh
# After zipping, removing the redundant text files  
# 2019-03-29

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/merged_softlinks_textfile"

rm  $inmain/*.txt
