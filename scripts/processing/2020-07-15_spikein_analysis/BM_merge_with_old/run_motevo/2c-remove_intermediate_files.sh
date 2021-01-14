#!/bin/sh
# Jake Yeung
# 2c-remove_intermediate_files.sh
#  
# 2020-11-09

jmark="H3K27me3"
# jmark="H3K4me3"
maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged3_Peaks"

indir=${maindir}/${jmark}
# cd ${maindir}/${jmark}

fdir=${indir}/fasta
fdirsplit=${indir}/fastasplit

[[ ! -d $fdir ]] && echo "$fdir not found, exiting" && exit 1
[[ ! -d $fdirsplit ]] && echo "$fdirsplit not found, exiting" && exit 1

# check tar is available
[[ ! -e $indir/fasta.tar.gz ]] && echo "$indir/fasta.tar.gz not found, exiting" && exit 1
[[ ! -e $indir/fastasplit.tar.gz ]] && echo "$indir/fastasplit.tar.gz not found, exiting" && exit 1

# gigive command rto remove dir
cmd="rm -r $fdir $fdirsplit"
echo $cmd

# remove motevo_output dir
indir2=${indir}/motevo_outputs
[[ ! -d $indir2 ]] && echo "$indir2 not found, exiting" && exit 1
cd $indir2

d1=${indir2}/split
d2=${indir2}/merged
d3=${indir2}/closestbed_multiple_genes

# check tars

[[ ! -e $indir2/split.tar.gz ]] && echo "$indir2/split.tar.gz not found, exiting" && exit 1
[[ ! -e $indir2/merged.tar.gz ]] && echo "$indir2/merged.tar.gz not found, exiting" && exit 1
[[ ! -e $indir2/closestbed_multiple_genes.tar.gz ]] && echo "$indir2/closestbed_multiple_genes.tar.gz not found, exiting" && exit 1

cmd="rm -r $d1 $d2 $d3"
echo $cmd
