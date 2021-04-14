#!/bin/sh
# Jake Yeung
# 2-copy_H3K27me3_fastq_to_hpc.sh
#  
# 2020-11-21

indir="/hpc/archive/hub_oudenaarden/seqdata/VAN5046/200910_NS500813_0646_AHFLJGBGXG/Data/Intensities/BaseCalls/AVOEI856-27"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged"

cd $indir

for f in $(find . -name "PZ-ChIC-mouse-BM-H3K27me3*fastq.gz" 2>/dev/null); do
    cp $f $outdir
done
