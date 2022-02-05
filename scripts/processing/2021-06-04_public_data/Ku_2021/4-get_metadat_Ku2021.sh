#!/bin/sh
# Jake Yeung
# 4-get_metadata_Bartosovic.sh
#  
# 2021-06-07


# wait 10000

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data"
# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/metadata_Bartosovic.txt"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/metadata_Ku2021.txt"
# db="/hpc/hub_oudenaarden/jyeung/data/public_data/SRAmetadb.sqlite"
# pysradb metadata SRP227999 --saveto $outf --db $db
# pysradb metadata GSE139857 --saveto $outf --db $db
pysradb metadata SRP227999 --saveto $outf


