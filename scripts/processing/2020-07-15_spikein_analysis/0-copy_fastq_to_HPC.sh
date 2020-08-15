#!/bin/sh
# Jake Yeung
# 1-copy_fastq_to_HPC.sh
#  
# 2020-07-15

indir="/hpc/hub_oudenaarden/seqdata/VAN4786/200714_NS500813_0627_AH25KMBGXG/Data/Intensities/BaseCalls/AVOEI856-8"
cd $indir

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs"

find . -name "*PZ*.fastq*" -exec cp {} $outdir \; 
