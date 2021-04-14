#!/bin/sh
# Jake Yeung
# 1-run_demux.sh
# 2020-08-10

jmem='16G'
jtime='4:00:00'

# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5230_BM"
# inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM/mouse"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
BNAME=$inmain/demux
[[ ! -d $BNAME ]] && mkdir $BNAME
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters  # important
demux.py *.fastq.gz -hd 0 | grep submission | sh
