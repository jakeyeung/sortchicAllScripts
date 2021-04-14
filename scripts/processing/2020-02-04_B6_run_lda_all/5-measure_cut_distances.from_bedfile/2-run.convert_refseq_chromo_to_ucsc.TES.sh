#!/bin/sh
# Jake Yeung
# convert_refseq_chromo_to_ucsc.sh
#  
# 2020-08-06

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/5-measure_cut_distances.from_bedfile/convert_refseq_chromo_to_ucsc.py"

# inf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.sga.gz"
inf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTes.sga.gz"
reffile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GCF_000001635.26_GRCm38.p6_assembly_report.txt"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTes.sga.chromorenamed.gz"

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -assembly_report $reffile -outfile $outf --add_chr
