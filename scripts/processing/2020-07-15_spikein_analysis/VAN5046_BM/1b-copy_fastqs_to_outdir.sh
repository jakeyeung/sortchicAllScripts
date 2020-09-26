#!/bin/sh
# Jake Yeung
# 1b-copy_fastqs_to_outdir.sh
#  
# 2020-08-10

# outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN6969"
# indir="/hpc/hub_oudenaarden/seqdata/VAN4969/200807_NS500813_0635_AHTWVNBGXF/Data/Intensities/BaseCalls/AVOEI856-15"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/200904_NS500813_0644_AHFLC3BGXG/Data/Intensities/BaseCalls/AVOEI856-20"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/200910_NS500813_0646_AHFLJGBGXG/Data/Intensities/BaseCalls/AVOEI856-27"

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM"
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

cd $indir
find . -name "*PZ*.fastq*" -exec cp {} $outdir \;
