#!/bin/sh
# Jake Yeung
# create_genomic_windows.sh
# Make genomic windows
# 2018-12-14
# From http://quinlanlab.org/tutorials/bedtools/answers.html

winsize=5000
stepsize=1000
chromsizes="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/chromsizes.mm10.txt"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/genomebins/mm10.winsize.$winsize.stepsize.$stepsize.bed"

bedtools makewindows -g $chromsizes -w 5000 -s 1000 > $outf
