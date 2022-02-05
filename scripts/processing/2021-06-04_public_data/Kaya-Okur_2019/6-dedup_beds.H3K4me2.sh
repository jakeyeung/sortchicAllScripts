#!/bin/sh
# Jake Yeung
# 6-dedup_beds.sh
#  
# 2021-06-07

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K4me2"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K4me2/dedup_beds"
[[ ! -d $outdir ]] && mkdir $outdir

for b in `ls -d ${indir}/*.bed.gz`; do
    bmain=$(basename $b)
    bmain=${bmain%%.*}
    bout=$outdir/${bmain}.dedup.bed
    # echo $bout
    zcat $b | awk '!a[$0]++' > $bout
    # exit 0
done
