#!/bin/sh
# Jake Yeung
# copy_blacklist.sh
#  
# 2019-12-03

inf="/home/jyeung/hpc/databases/blacklists/mm10.blacklist.copy.bed.gz"
# outf="~/data/scchic/databases/mm10.blacklist.bed.gz"
outf="/home/jyeung/data/scchic/databases/blacklists/mm10.blacklist.bed.gz"

cp $inf $outf
