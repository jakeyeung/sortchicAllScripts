#!/bin/sh
# Jake Yeung
# 1-run_demux.sh
# 2021-10-29

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs"

jmem='16'  # no suffix for SCMO
jtime='24'
BNAME=$inmain/demux
[[ ! -d $BNAME ]] && mkdir $BNAME
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate SCMO_2021
demux.py *.fastq.gz -hd 0 --cluster -mem $jmem -time $jtime -sched slurm
