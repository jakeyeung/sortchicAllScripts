#!/bin/sh
# Jake Yeung
# 3b-bam_to_bigwig_with_scale_factors.sh
#  
# 2020-06-19

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/2-make_bigwigs/FewerClusters/calculate_scale_factors_for_deeptools.py"

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits${suffix}"

inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits.Downsamp/H3K4me3-BM_AllMerged.HSPCs.sorted.cleaned.bam"

# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inbam

python $ps -infile $inbam
