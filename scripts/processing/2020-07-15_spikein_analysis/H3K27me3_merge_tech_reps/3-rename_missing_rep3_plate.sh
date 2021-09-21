#!/bin/sh
# Jake Yeung
# 3-rename_missing_rep3_plate.sh
# plates 10, 4, 6 in rep3 have different names. Rename them
# 2020-11-21

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged"

cd $indir
# PZ-ChIC-mouse-BM-H3K27me3-10_AHFLJGBGXG_S3_L002_R2_001.fastq.gz
for f in `ls -d PZ-ChIC-mouse-BM-H3K27me3-10_*.fastq.gz`; do
    p1=$(echo $f | cut -d"_" -f1)
    p2=$(echo $f | cut --complement -d"_" -f1)
    # replace p1 with PZ-BM-rep3-H3K27me3-10
    p1new="PZ-BM-rep3-H3K27me3-10_"
    fnew=${p1new}${p2}
    # echo $f ">" $fnew
    mv $f $fnew
done

# PZ-ChIC-mouse-BM-H3K27me3-6_AHFLJGBGXG_S2_L003_R2_001.fastq.gz 
for f in `ls -d PZ-ChIC-mouse-BM-H3K27me3-6_*.fastq.gz`; do
    p1=$(echo $f | cut -d"_" -f1)
    p2=$(echo $f | cut --complement -d"_" -f1)
    # replace p1 with PZ-BM-rep3-H3K27me3-6
    p1new="PZ-BM-rep3-H3K27me3-6_"
    fnew=${p1new}${p2}
    # echo $f ">" $fnew
    mv $f $fnew
done

# PZ-ChIC-mouse-BM-H3K27me3-4_AHFLJGBGXG_S1_L004_R2_001.fastq.gz
for f in `ls -d PZ-ChIC-mouse-BM-H3K27me3-4_*.fastq.gz`; do
    p1=$(echo $f | cut -d"_" -f1)
    p2=$(echo $f | cut --complement -d"_" -f1)
    # replace p1 with PZ-BM-rep3-H3K27me3-4
    p1new="PZ-BM-rep3-H3K27me3-4_"
    fnew=${p1new}${p2}
    # echo $f ">" $fnew
    mv $f $fnew
done
