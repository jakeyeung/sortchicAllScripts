#!/bin/sh
# Jake Yeung
# copy_blacklist.sh
#  
# 2022-01-24

inf="t2:/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.nospikeins.nochromo.bed"

outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases"

scp $inf $outdir
