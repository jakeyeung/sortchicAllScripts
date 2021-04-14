#!/bin/sh
# Jake Yeung
# 0-download_refseq_annots.sh
#  
# 2020-08-06

# link="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.gff.gz"
# link="GCF_000002035.6_GRCz11_genomic.gbff.gz"
link="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.gbff.gz"
outdir="/hpc/hub_oudenaarden/jyeung/data/databases/zf/refseq"
cd $outdir
curl -O $link
