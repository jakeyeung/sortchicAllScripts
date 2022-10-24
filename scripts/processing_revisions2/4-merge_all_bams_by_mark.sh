#!/bin/sh
# Jake Yeung
# 4-merge_all_bams_by_mark.sh
#  
# 2022-07-19

jmem='96G'
jtime='24:00:00'

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged"

outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/old_new_bams_merged"

module load samtools
for mark in $marks; do
		echo $mark
		outf=${outdir}/BM_allmerged_${mark}.bam
		infs=${inmain}/${mark}/*.bam
		echo $infs
		cmd="module load samtools; samtools merge $outf $infs; samtools index $outf"

		BNAME=$outdir/sbatch_output_${mark}
		DBASE=$(dirname "${BNAME}")
		[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

		sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=merge_${mark} --wrap "$cmd"
done
