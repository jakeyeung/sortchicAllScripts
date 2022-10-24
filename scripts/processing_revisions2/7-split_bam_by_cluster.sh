#!/bin/sh
# Jake Yeung
# run.split_bam_by_cluster.sh
#  
# 2020-01-09

jmem='32G'
jtime='24:00:00'

ps="/nfs/scistore12/hpcgrp/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
mapq=40

indir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/old_new_bams_merged"

annotdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/metadata"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $annotdir ]] && echo "$annotdir not found, exiting" && exit 1

outmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/old_new_bams_merged/split_by_celltype"

queue="gpu"
for jmark in $jmarks; do
	outdir=${outmain}/${jmark}
	[[ ! -d $outdir ]] && mkdir $outdir
	[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1  # make dir before hand
	bname="BM_allmerged_${jmark}"
    inf=${indir}/${bname}.bam
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

	annotname="cell_ctype_metadata_for_splitting_bams.${jmark}.txt"
    annotfile=${annotdir}/${annotname}
    [[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
	module load R/4.1.2
	cmd=". /nfs/scistore12/hpcgrp/jyeung/miniconda3/etc/profile.d/conda.sh; conda activate scmo2022; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq ${mapq} --add_chr_prefix --annot_no_colnames"
	sbatch -p ${queue} --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=split_${jmark} --wrap "$cmd"
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo2022; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq ${mapq} --add_chr_prefix" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N split_${jmark} -m beas -M j.yeung@hubrecht.eu
done

