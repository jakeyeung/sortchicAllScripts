#!/bin/sh
# Jake Yeung
# 4-bam_to_bigwig.sh
#  
# 2019-02-04


. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

bamCoverage
