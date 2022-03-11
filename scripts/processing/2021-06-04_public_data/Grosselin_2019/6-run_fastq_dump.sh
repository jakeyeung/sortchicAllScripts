#!/bin/sh
# Jake Yeung
# 6-run_fastq_dump.sh
#  
# 2021-06-20

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Grosselin_et_al_2019/SRA_data/prefetch_outputs/SRR7536859"

cd $indir

fastq-dump SRR7536859.sra
