#!/bin/sh
# Jake Yeung
# 3-rename_fastq_files.sh
#  
# 2020-11-21

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged"

# part1 (first seq) are named: PZ-BM-rep3-H3K9me3-H3K4me1-3_HN7YKBGXG_S15_L004_R1_001.fastq.gz
# part2 (reseq) are named: PZ-ChIC-BM-rep3-H3K9me3-H3K4me1-1_AH22Y7BGXH_S15_L004_R1_001.fastq.gz

# rename reseq to match firstseq
cd $indir
echo $indir

for f in `ls -d PZ-ChIC-BM-rep3-H3K9me3-H3K4me1-*fastq.gz`; do
    # PZ-ChIC-BM-rep3-H3K9me3-H3K4me1-4_AH22Y7BGXH_S18_L003_R1_001.fastq.gz
    # echo $f
    p1=$(echo $f | cut -d"-" -f1-6)
    p2=$(echo $f | cut --complement -d"-" -f1-6)
    # echo "${p1}-${p2}" "==" $f
    p1new="PZ-BM-rep3-H3K9me3-H3K4me1"
    fnew="${p1new}-${p2}"
    echo $f "->" $fnew
    mv $f $fnew
done
