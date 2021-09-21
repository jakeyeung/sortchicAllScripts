#!/bin/sh
# Jake Yeung
# 0-download_assembly_reports.sh
#  
# 2020-08-06

link="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Danio_rerio/latest_assembly_versions/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_assembly_report.txt"

outdir="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq"

cd $outdir

curl -O $link
