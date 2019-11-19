#!/bin/sh
# Jake Yeung
# 4-run_fastqc.sh
#  
# 2019-09-28

n=0
maxjobs=4

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3985/raw_demultiplexed"

jmem='4G'
jtime='12:00:00'

for indir in `ls -d $inmain/PZ*`; do 
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
    outdir=$indir/fastqc
    [[ ! -d $outdir ]] && mkdir $outdir
    fastqc $indir/*.fastq.gz -o $outdir&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
