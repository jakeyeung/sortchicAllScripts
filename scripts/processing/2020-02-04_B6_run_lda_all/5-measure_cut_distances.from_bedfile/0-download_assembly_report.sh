#!/bin/sh
# Jake Yeung
# 0-download_assembly_report.sh
#  
# 2020-08-06

link="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_assembly_report.txt"

outdir="/hpc/hub_oudenaarden/jyeung/data/databases/refseq"
cd $outdir

curl -O $link
