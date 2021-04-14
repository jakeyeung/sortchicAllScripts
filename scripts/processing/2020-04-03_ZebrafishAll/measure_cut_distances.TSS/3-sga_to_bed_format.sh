#!/bin/sh
# Jake Yeung
# 3-sga_to_bed_format.sh
#  
# 2020-08-06

# inf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.sga.chromorenamed.gz"
inf="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.parsed.parsed.TSS.chromorenamed.txt.gz"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq/GCF_000002035.6_GRCz11_genomic.parsed.parsed.TSS.chromorenamed.bed.gz"

zcat $inf | awk 'BEGIN{OFS="\t"}{print ($1, $3, ($3 + 1), $4, $6)}' | gzip > $outf
