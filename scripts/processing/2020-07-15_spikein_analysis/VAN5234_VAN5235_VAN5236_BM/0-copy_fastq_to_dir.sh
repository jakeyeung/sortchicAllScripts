#!/bin/sh
# Jake Yeung
# 0-copy_fastq_to_dir.sh
#  
# 2020-10-03

# inmain="/hpc/hub_oudenaarden/nborst/201001_NS500414_0731_HNHHHBGXG/Data/Intensities/BaseCalls/AVOU856-40"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM/201002_NS500813_0658_AHNGVJBGXG/Data/Intensities/BaseCalls/AVOU856-42"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5230_BM"

cd $inmain

find  . -name "PZ*.fastq.gz" -exec cp {} ${outdir} \; 
