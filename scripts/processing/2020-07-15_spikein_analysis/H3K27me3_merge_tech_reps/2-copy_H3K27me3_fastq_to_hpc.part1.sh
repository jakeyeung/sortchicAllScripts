#!/bin/sh
# Jake Yeung
# 2-copy_H3K27me3_fastq_to_hpc.sh
# copies all fastq except for rep10, 4, 6 which are found in VAN5046 
# without rep3 names so better to copy them separately
#  
# 2020-11-21

indir="/hpc/archive/hub_oudenaarden/seqdata"

fastqlist="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/H3K27me3_merge_tech_reps/list_H3K27me3_reps_to_merge.again.out"  # relative to indir (/hpc/archive/hub_oudenaarden/seqdata)
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged"

# read fastq
while read fname; do
  # echo "${fname}"
  inf=${indir}/${fname}
  cp $inf $outdir
done < $fastqlist
