#!/bin/sh
# Jake Yeung
# 3-link_old_new_bams_to_one_dir.sh
#  
# 2022-07-19

declare -A markray

markray[H3K4me1]=k4me1
markray[H3K4me3]=k4me3
markray[H3K27me3]=k27me3
markray[H3K9me3]=k9me3

outmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged"
# inmainnew="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/tagged_bams"
# inmainold="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/bams"

for mark in "${!markray[@]}"; do
    printf "%s is in %s\n" "$mark" "${markray[$mark]}"
	outdir=${outmain}/$mark
	cd $outdir

    inmainnew="../../new_experiments/raw_data/tagged_bams"
    inmainold="../../first_submission_data/bams"

	indirnew="${inmainnew}/BM_${markray[$mark]}"
	indirold=${inmainold}/$mark
	[[ ! -d $indirnew ]] && echo "$indirnew not found, exiting" && exit 1
	[[ ! -d $indirold ]] && echo "$indirold not found, exiting" && exit 1

	echo "fnews"
	for fnew in `ls $indirnew/*.bam*`; do 
			echo $fnew
			ln -s $fnew $outdir
	done

	echo "folds"
	for fold in `ls $indirold/*.bam*`; do 
			echo $fold
			ln -s $fold $outdir
	done

done

# 
# for mark in $marks; do
# 
# 	
# 
# done
