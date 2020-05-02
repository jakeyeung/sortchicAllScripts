#!/bin/sh
# Jake Yeung
# 2-make_blacklist_with_filtered_bins.sh
#  
# 2020-01-10

# mapq=40
# [[ $mapq != [0-9]* ]] && echo "Must be integer: $mapq" && exit 1
#  --ignoreDuplicates --minMappingQuality $mapq"
# 
bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.bed.gz"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1
correlatedbins="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/from_rstudioiserver/quality_controls_intestines.OK/correlated_bins.Scraped.bincutoff_0.95.2019-12-22.bed"
[[ ! -e $correlatedbins ]] && echo "$correlatedbins not found, exiting" && exit 1

qcdir="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/from_rstudioiserver/quality_controls_intestines.OK"
outf="${qcdir}/blacklist_corr_bins_merged.bed"
[[ -e $outf ]] && echo "$outf found, exiting" && exit 1
outfmerged="${qcdir}/blacklist_corr_bins_merged.merged_sorted.bed"

zcat $bl | sort -k1,1 -k2,2n >> $outf
# remove chrchr -> chr
sed 's/chrchr/chr/g' $correlatedbins | sort -k1,1 -k2,2n >> $outf

sort --output $outf -k1,1 -k2,2n $outf
bedtools merge -i $outf > ${outfmerged}

ret=$?; [[ $ret -ne 0  ]] && echo "ERROR: script failed" && exit 1

# clean up 
rm $outf

echo "Good blacklist is here:"
echo ${outfmerged}
