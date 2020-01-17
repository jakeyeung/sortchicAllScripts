#!/bin/sh
# Jake Yeung
# demultiplex_fastqs.sh
#  
# 2019-09-28

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6/AVO479"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6.demuxbugfixed"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6.quotebugfixed"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6.quotebugfixed"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6"
grepstr="*.fastq.gz"

jmem='8G'
jtime='2:00:00'
BNAME=$inmain/demux
[[ ! -d $BNAME ]] && mkdir $BNAME
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

# demux.py *.fastq.gz -hd 0 | grep submission | sh
demux.py *.fastq.gz -hd 0 -merge _ | grep submission | sh

