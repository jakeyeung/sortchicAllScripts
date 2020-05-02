#!/bin/sh
# Jake Yeung
# 0-make_links_from_archive_all.sh
#  
# 2019-12-13

# copy everything

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_linneg_K4me1"
[[ ! -d $outdir ]] && mkdir $outdir

indir="/hpc/archive/hub_oudenaarden/userBackups/from_sequencer/200120_NS500413_0666_AH2GM5BGXF/Data/Intensities/BaseCalls"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
cd $indir

find . -name "*PZ-ChIC-Bl6-BM-lin-H3K4me1*.fastq.gz" 2>/dev/null | xargs -I {} cp {} $outdir/.  # the original wildtype data
