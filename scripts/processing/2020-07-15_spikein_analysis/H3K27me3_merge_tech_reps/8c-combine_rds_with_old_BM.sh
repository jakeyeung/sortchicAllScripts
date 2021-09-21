#!/bin/sh
# Jake Yeung
# 8c-combine_rds_with_old_BM.sh
#  
# 2020-11-28

jmem='32G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/H3K27me3_merge_tech_reps/merge_rds_mats.R"

inf1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables/for_LDA/PZ-BM-rep2and3-H3K27me3-rep3Reseq.binsize_50000.rds"
inf2="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs/count_mat_cbfilt_maxcountsfilt.all.2020-06-16.H3K27me3.rds"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables/for_LDA"
outf="${outdir}/PZ-BM-rep2and3withold-H3K27me3-rep3Reseq.binsize_50000.rds"

BNAME=${outdir}/combine_new_old_reseeq
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $inf1 $inf2 -outfile $outf"

sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=combine_new_old_reseq --wrap "$cmd"
