#!/bin/sh
# Jake Yeung
# 3-run_mat_to_mat_long.sh
# Make format so you easily put into sql databsae (i.e., long format)
# 2016-08-03

longscript="/Home/jyeung/projects/sleep_deprivation/atacseq_analysis/2-extract_dhs_signal/mat_to_long_mat.R"
# inf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/covbed.atacseq/dhs_signal_windows500.long.mat"
# outf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/covbed.atacseq/dhs_signal_windows500.long.longforsql.mat"
inf="/archive/epfl/upnae/jyeung/sleep_deprivation/atacseq_signal/covbed.unmerged.atacseq/atacseq_signal_windows500.long.mat"
outf="/archive/epfl/upnae/jyeung/sleep_deprivation/atacseq_signal/covbed.unmerged.atacseq/atacseq_signal_windows500.long.longforsql.mat"

# metadat="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/covbed.atacseq/metadata.txt"
metadat="/archive/epfl/upnae/jyeung/sleep_deprivation/atacseq_signal/covbed.unmerged.atacseq/metadata.txt"  # do samples really match bed files?
cnames_samps=`cat $metadat | grep -Po '(?<=atacseq/).+?(?=.bam)' | tr '\n' ',' | sed 's/,$//g'`
cnames="chromo,start,end,gene,dist,"$cnames_samps
# replace characters between commas to integer
# echo $cnames_samps
# echo $metadata
nsamps=$(cat $metadat | wc -l)
# nsamps="cat $metadat | awk 'BEGIN{FS=","} {print NF?NF-1:0}"
classes_samps=`for f in $(seq $nsamps); do echo "integer"; done | tr "\n" "," | sed 's/.$//'`
# for f in $(seq 40); do echo "integer"; done | wc -l
# classes_samps=`echo $cnames_samps | sed 's/[^,]/X/g' | sed 's/XXXXXXXX\|XXXXXXX/integer/g'`
# echo "Classes: $classes_samps"
# 															  XXXXXXX
# classes_samps=`echo $cnames_samps | sed 's/[^,]/X/g'`
# classes_samps=`echo $cnames_samps | tr -s '[^,]' 'X'`
classes="character,integer,integer,character,integer,"$classes_samps


[[ ! -e $longscript ]] && echo "$longscript not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

# echo $cnames_samps
# echo $cnames
# echo $classes

echo "Rscript $longscript $inf $outf $cnames $classes"
Rscript $longscript $inf $outf $cnames $classes
