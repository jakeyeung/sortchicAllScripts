#!/bin/sh
# Jake Yeung
# 2-copy_H3K27me3_fastq_to_hpc.sh
#  
# 2020-11-21

# indir="/hpc/archive/hub_oudenaarden/seqdata/VAN5046/200910_NS500813_0646_AHFLJGBGXG/Data/Intensities/BaseCalls/AVOEI856-27"
indir="/hpc/hub_oudenaarden/Peter/data/VAN5277"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged"

cd $indir

for f in $(find . -name "PZ-ChIC-BM-rep3-H3K9me3-H3K4me1*fastq.gz" 2>/dev/null); do
    cp $f $outdir
done
