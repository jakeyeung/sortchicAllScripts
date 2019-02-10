#!/bin/sh
# Jake Yeung
# 4-cleanup_files.sh
# Clean up files after generating what you need 
# 2019-02-03

# dirmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output"
dirmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_multigene/H3K4me1"

# dnames="fasta fastasplit tmp motevo_outputs/split motevo_outputs/merged motevo_outputs/closestbed_multiple_genes"
dnames="fasta fastasplit motevo_outputs/split motevo_outputs/merged motevo_outputs/closestbed_multiple_genes"

for d in $dnames; do
# for d in $(ls -d $dirmain/$dnames); do
    jdir=$dirmain/$d
    [[ ! -d $jdir ]] && echo "$jdir not found, exiting" && exit 1
    # rm -r  $jdir
done


