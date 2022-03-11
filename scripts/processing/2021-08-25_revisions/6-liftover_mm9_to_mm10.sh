#!/bin/sh
# Jake Yeung
# 6-liftover_mm9_to_mm10.sh
#  
# 2021-08-27

liftOver="/hpc/hub_oudenaarden/jyeung/software/ucsc_utils/liftOver"
mapchain="/hpc/hub_oudenaarden/jyeung/data/databases/chainfiles/mm9ToMm10.over.chain.gz"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.bonemarrow_filt.rownames_filt.mm9.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.bonemarrow_filt.rownames_filt.mm10.bed"
unmapped="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data/atac_matrix.bonemarrow_filt.rownames_filt.mm10.unmapped"

$liftOver $inf $mapchain $outf $unmapped
