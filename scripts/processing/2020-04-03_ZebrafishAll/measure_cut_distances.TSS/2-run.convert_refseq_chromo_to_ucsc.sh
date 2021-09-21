#!/bin/sh
# Jake Yeung
# convert_refseq_chromo_to_ucsc.sh
#  
# 2020-08-06

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/5-measure_cut_distances.from_bedfile/convert_refseq_chromo_to_ucsc.py"

# inf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.sga.gz"
# reffile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GCF_000001635.26_GRCm38.p6_assembly_report.txt"
# outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.sga.chromorenamed.gz"
inf="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.parsed.parsed.TSS.txt.gz"
reffile="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_assembly_report.txt"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.parsed.parsed.TSS.chromorenamed.txt.gz"

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -assembly_report $reffile -outfile $outf --add_chr
