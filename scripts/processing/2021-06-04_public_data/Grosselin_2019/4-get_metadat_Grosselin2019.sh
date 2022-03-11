#!/bin/sh
# Jake Yeung
# 4-get_metadata_Bartosovic.sh
#  
# 2021-06-07


# wait 10000

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Grosselin_et_al_2019/SRA_data"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Grosselin_et_al_2019/SRA_data/metadata_Grosselin2019.txt"
srp="SRP154341"
cd $indir
pysradb metadata $srp --saveto $outf


