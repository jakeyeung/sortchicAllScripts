#!/bin/sh
# Jake Yeung
# 1b-copy_fastqs_to_outdir.sh
#  
# 2020-08-10

# outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN6969"
indir="/hpc/hub_oudenaarden/seqdata/VAN4969/200807_NS500813_0635_AHTWVNBGXF/Data/Intensities/BaseCalls/AVOEI856-15"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN6969/K562"

cd $indir

find . -name "*K562*.fastq*" -exec cp {} $outdir \;
