#!/bin/sh
# Jake Yeung
# 6c-expand_multigene_to_each_line_dhssignal.sh
# Expand @ delim genes to individual rows
# 2016-04-05

decompress="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/merged_dhs/convert_compressed_bed_to_long.py"
# inbed="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/covbed.atacseq/dhs_signal_windows500.mat"  # merged bed
# outbed="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/covbed.atacseq/dhs_signal_windows500.long.mat"
inbed="/archive/epfl/upnae/jyeung/sleep_deprivation/atacseq_signal/covbed.unmerged.atacseq/atacseq_signal_windows500.mat"  # unmerged bed
outbed="/archive/epfl/upnae/jyeung/sleep_deprivation/atacseq_signal/covbed.unmerged.atacseq/atacseq_signal_windows500.long.mat"  # unmerged bed

[[ ! -e $decompress ]] && echo "$decompress not found, exiting" && exit 1
[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1

python $decompress $inbed $outbed --has_dist --gene_col_i 3
