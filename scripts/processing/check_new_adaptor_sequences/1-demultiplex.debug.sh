#!/bin/sh
# Jake Yeung
# 1-demultiplex.sh
#  
# 2019-09-04

# inmain="/hpc/hub_oudenaarden/mvanins/newruns/oud3697/190830_NS500414_0628_H7YVNBGXC/Data/Intensities/BaseCalls/AVO505"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/from_vanins/AVO505"
# inmain="/hpc/hub_oudenaarden/bdebarbanson/data/jake"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/fastq.debugging"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3697/raw_demultiplexed.fastqOneDir2"
# [[ ! -d $outdir ]] && mkdir $outdir

# grepstr="OUD3697A*/PZ-ChIC*.fastq.gz"
grepstr="PZ-ChIC*.fastq.gz"

# cd $inmain

# jmem='32G'
# jtime='24:00:00'
# BNAME=$outdir/demux.qsub
# DBASE=$(dirname "${BNAME}")
# [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1


cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
echo "demux.py "${grepstr}" -o $outdir -merge _ | grep submission | sh"
