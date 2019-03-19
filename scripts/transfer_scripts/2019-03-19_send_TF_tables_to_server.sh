#!/bin/sh
# Jake Yeung
# 2019-03-19_send_TF_tables_to_server.sh
#  
# 2019-03-19

inf="/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects.RData"
outf="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/."

scp $inf $outf
