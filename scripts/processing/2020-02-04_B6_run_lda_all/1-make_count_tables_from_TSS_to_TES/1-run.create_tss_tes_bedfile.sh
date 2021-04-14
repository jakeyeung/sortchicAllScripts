#!/bin/sh
# Jake Yeung
# 1-run.create_tss_tes_bedfile.sh
#  
# 2020-08-10

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/1-make_count_tables_from_TSS_to_TES/create_tss_tes_bedfile.py"

jstrand="neg"
inf1="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.${jstrand}.bed.gz"
inf2="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTes.chromorenamed.${jstrand}.bed.gz"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTssToTes.chromorenamed.${jstrand}.bed.gz"

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -TssBed $inf1 -TesBed $inf2 -OutBed $outf

