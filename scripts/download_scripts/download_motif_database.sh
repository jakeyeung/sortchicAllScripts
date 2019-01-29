#!/bin/sh
# Jake Yeung
# download_motif_database.sh
# Download big motif database 
# 2019-01-14

outdir="/hpc/hub_oudenaarden/jyeung/data/databases/motifs"
cd $outdir

# Specify database name:
# feather_database_url='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc8nr/region_based/hg19-regions-9species.all_regions.mc8nr.feather'
# feather_database_url='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/region_based/hg19-regions-1M-9species.all_regions.mc9nr.feather'
feather_database_url='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc8nr/region_based/hg19-regions-9species.all_regions.mc8nr.feather'

feather_database="${feather_database_url##*/}"

# Download database directly (with wget or curl):
# wget "${feather_database_url}"
curl -O "${feather_database_url}"

# Download sha256sum.txt (with wget or curl):
# wget https://resources.aertslab.org/cistarget/databases/sha256sum.txt
# curl -O https://resources.aertslab.org/cistarget/databases/sha256sum.txt

# Check if sha256 checksum matches for the downloaded database:
# awk -v feather_database=${feather_database} '$2 == feather_database' sha256sum.txt | sha256sum -c -

# If you downloaded mulitple databases, you can check them all at onces with:
# sha256sum -c sha256sum.txt
