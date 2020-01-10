#!/bin/sh
# Jake Yeung
# 5-make_links_to_tables.sh
#  
# 2019-11-08

indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3910_HD_0/countTables_geneTSS"
indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_HD_0/countTables_geneTSS"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_oud3910_merged/countTables_geneTSS"

ln -s $indir1/*.RData $outdir
ln -s $indir2/*.RData $outdir
