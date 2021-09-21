#!/bin/sh
# Jake Yeung
# 2-copy_H3K9me3-H3K4me1_part1.sh 
# copy first seequencing run of double stain. Second run not yet uploaded in seqdata, it is still in Peters folder
# 2020-11-21

indir="/hpc/archive/hub_oudenaarden/seqdata"

# fastqlist="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/H3K27me3_merge_tech_reps/list_H3K27me3_reps_to_merge.again.out"  # relative to indir (/hpc/archive/hub_oudenaarden/seqdata)
# fastqlist="list_H3K4me1-H3K9me3_reps_to_merge.again.out"

fastqlist="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/H3K4me1-H3K9me3_dbl_merge_tech_reps/list_H3K4me1-H3K9me3_reps_to_merge.again.out"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged"

# read fastq
while read fname; do
  # echo "${fname}"
  inf=${indir}/${fname}
  cp $inf $outdir
done < $fastqlist
