#!/bin/sh
# Jake Yeung
# 1-demultiplex.sh
#  
# 2019-09-04

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/fastq"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3697/raw_demultiplexed.nomerge"

grepstr="OUD3697A*/PZ-ChIC*.fastq.gz"

# cd $inmain

# jmem='32G'
# jtime='24:00:00'
# BNAME=$outdir/demux.qsub
# DBASE=$(dirname "${BNAME}")
# [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1


cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

# demux.py "$inmain/${grepstr}" -o $outdir -hd 1 -merge _ | grep submi | sh 
# demux.py "$inmain/${grepstr}" -o $outdir | grep submission | sh 
demux.py "$inmain/${grepstr}" -o $outdir | grep submission | sh 
# demux.py "$inmain/${grepstr}" -o $outdir -merge _ | grep submission | sh 
