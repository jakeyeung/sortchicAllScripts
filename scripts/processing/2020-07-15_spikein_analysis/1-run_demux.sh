#!/bin/sh
# Jake Yeung
# 1-run_demux.sh
#  
# 2019-12-19

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataK562"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs"

jmem='16G'
jtime='4:00:00'
BNAME=$inmain/demux
[[ ! -d $BNAME ]] && mkdir $BNAME
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo3
demux.py *.fastq.gz -hd 0 | grep submission | sh
# demux.py *.fastq.gz -hd 0
