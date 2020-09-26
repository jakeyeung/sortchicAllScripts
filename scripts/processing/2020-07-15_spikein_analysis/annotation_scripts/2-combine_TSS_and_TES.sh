#!/bin/sh
# Jake Yeung
# 2-combine_TSS_and_TES.sh
#  
# 2020-08-14

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/1-make_count_tables_from_TSS_to_TES/create_tss_tes_bedfile.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

# jstrands="pos neg"

# for jstrand in $jstrands; do
    inf1="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS.txt.bed.gz"
    inf2="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TES.txt.bed.gz"
    outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS_TES.txt.bed.gz"
    python $ps -TssBed $inf1 -TesBed $inf2 -OutBed $outf 
# done
