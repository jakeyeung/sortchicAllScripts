#!/bin/sh
# Jake Yeung
# 6e-check_each_chromos_have_peaks.sh
#  
# 2020-11-02

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_60_maxcount_40"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_60_maxcount_40"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/merged_bams.first_and_second_rounds/hiddendomains_outputs_minlength_2500.FromR.maxcount_40_60_80"

echo $inmain

nchromoscheck="21"
for indir in `ls -d $inmain/*cutoff`; do
    bname=$(basename $indir)
    dname=$(echo $bname | cut --complement -d"." -f1)
    inf=${indir}/${dname}_analysis.bed
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    nchromos=`cut -f1 $inf | sort | uniq | wc -l`
    [[ ! $nchromos -eq $nchromoscheck ]] && echo "Warning: ${bname}/${dname}_analysis.bed $nchromos != $nchromoscheck" && continue
done
