#!/bin/sh
# Jake Yeung
# 2019-02_03_copy_db.sh
# Copy sql database
# 2019-02-03

inf="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/sql/motevo_merged.closest.long.sqlite3"
# outf="~/data/scchic/databases/motevo_merged.closest.long.sqlite3"
outf="/Users/yeung/data/scchic/databases/"
tmpdir="/Users/yeung/data/scchic/databases/tmpdir"

[[ ! -d $tmpdir ]] && mkdir $tmpdir

rsync -avrL --partial-dir=$tmpdir --copy-links $inf $outf 
