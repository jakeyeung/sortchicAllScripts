#!/bin/sh
# Jake Yeung
# demultiplex_fastqs.sh
#  
# 2019-09-28

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3910_HD_0"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3910_HD_0"
# inmain="/hpc/hub_oudenaarden/fsalmen/fastq/191107_NS500413_0640_AHNM2HBGXC"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3985"

# do everything
# grepstr="PZ*.fastq.gz"
# do just oud3913 copied to oud3985
grepstr="PZ-ChIC-ZFWKM-H3K27me3-3_AHNJ35BGXC*.fastq.gz"

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

demux.py "${grepstr}" -hd 0 -merge _ | grep submission | sh

