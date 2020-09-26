#!/bin/sh
# Jake Yeung
# convert_refseq_chromo_to_ucsc.sh
#  
# 2020-08-06

jstrands="pos neg"

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/5-measure_cut_distances.from_bedfile/convert_refseq_chromo_to_ucsc.py"
reffile="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_assembly_report.txt"

    inf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS_TES.txt.bed.gz"
    outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS_TES.txt.chromorenamed.withchr.bed.gz"
    outf2="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS_TES.txt.chromorenamed.nochr.bed.gz"
    # outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS.txt.${jstrand}.chromorenamed.bed.gz"
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -assembly_report $reffile -outfile $outf --add_chr
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -assembly_report $reffile -outfile $outf2
