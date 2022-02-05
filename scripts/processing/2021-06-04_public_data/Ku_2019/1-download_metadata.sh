#!/bin/sh
# Jake Yeung
# 1-download_metadata.sh
#  
# 2021-06-09

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2019/SRA_data"
sra="SRP120016"

pysradb metadata $sra --saveto $outdir/metadata_Ku2019.txt
