#!/bin/sh
# Jake Yeung
# run.make_count_tables_filter_cells_for_LDA.sh
#  
# 2020-11-28

jmem='32G'
jtime='6:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-07-15_spikein_analysis/H3K27me3_merge_tech_reps/make_count_tables_filter_cells_for_LDA.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_tables"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outdir="${indir}/for_LDA"
[[ ! -d $outdir ]] && mkdir $outdir

annotfile="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.H3K27me3.2020-11-18.dupfilt.txt"
[[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1

cd $indir
infs=$(ls -d *50000.csv | tr '\n' ' ')
echo $infs
outf=${outdir}/"PZ-BM-rep3-H3K27me3-Round2Reseq.binsize_50000.rds"

BNAME=${outdir}/merge_filt_cells
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; 
cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -infile $infs -annotfile $annotfile -outfile $outf --add_chromo"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=merge_filt_cells --wrap "$cmd"
