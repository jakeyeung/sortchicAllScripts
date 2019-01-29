#!/bin/sh
# Jake Yeung
# 0-clean_fasta_database.sh
# Convert fasta database to uppercase to prevent MotEvo errors? 
# 2019-01-29

inf="/hpc/hub_oudenaarden/jyeung/data/databases/fasta/mm10.fa"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/fasta/mm10_upper.fa"

awk 'NR % 2 { print } !(NR % 2) {print toupper($0)}' $inf > $outf
