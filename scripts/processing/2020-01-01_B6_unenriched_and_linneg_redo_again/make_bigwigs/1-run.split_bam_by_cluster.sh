#!/bin/sh
# Jake Yeung
# run.split_bam_by_cluster.sh
#  
# 2020-01-09

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/H3K4me3-BM_Linneg_SC-merged.tagged.bam"
annotfile="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/from_rstudioiserver/pdfs_clustering/clusters_BM_All_merged_H3K4me3_20000_10000.txt"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_redo_2019-12-13.split_bams_by_clusters"

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq 40 --add_chr_prefix

