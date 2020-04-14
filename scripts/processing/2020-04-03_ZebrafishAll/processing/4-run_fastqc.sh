#!/bin/sh
# Jake Yeung
# 4-run_fastqc.sh
#  
# 2019-09-28


# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3910_HD_0/raw_demultiplexed"
inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/raw_data/raw_demultiplexed"

# jmem='4G'
# jtime='12:00:00'

for indir in `ls -d $inmain/PZ*`; do 
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
    outdir=$indir/fastqc
    [[ ! -d $outdir ]] && mkdir $outdir
    fastqc $indir/*.fastq.gz -o $outdir
done
