#!/bin/sh
# Jake Yeung
# 4-get_metadata_Bartosovic.sh
#  
# 2021-06-07


indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data"
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/metadata_Bartosovic.txt"
pysradb metadata SRP281199 --saveto $outf


