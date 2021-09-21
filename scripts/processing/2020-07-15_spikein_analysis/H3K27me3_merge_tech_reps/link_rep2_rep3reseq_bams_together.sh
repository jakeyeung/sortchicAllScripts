#!/bin/sh
# Jake Yeung
# link_rep2_rep3reseq_bams_together.sh
#  
# 2020-12-09

rep2dir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links"
rep3dir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/rep2_rep3reseq_bams_together"

ln -s $rep2dir/PZ-BM-rep2-H3K27me3*.bam* $outdir/.
ln -s $rep3dir/PZ-BM-rep3-H3K27me3*.bam* $outdir/.
