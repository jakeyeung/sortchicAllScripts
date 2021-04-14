#!/bin/sh
# Jake Yeung
# 12-split_bams_by_cluster.sh
#  
# 2020-12-11


jmem='8G'
jtime='1:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
jmarks="H3K27me3"
mapq=40

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BMround2all_VAN5046_VAN5109_VAN5230_BM_VAN5232_VAN5233_VAN5234_VAN5235_VAN5236/tagged_bams_links"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/rep2_rep3reseq_bams_together"
# annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables"
# annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios"
# annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04"
# annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged"
# annotfile="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.2020-12-10.txt"
# annotfile="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.2020-12-10.cluster_col_2.txt"
annotfile="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.2020-12-10.cluster_col_2.tab.txt"
[[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster"
[[ ! -d $outdir ]] && mkdir $outdir

[[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1

jmark="H3K27me3"
jmarkshort=`echo $jmark | sed 's/^H3//g'`  # in case of typoos H4K4me1

for inf in `ls -d ${indir}/PZ-BM-*${jmarkshort}*.bam`; do
    bname=$(basename $inf)
    bname=${bname%.*}

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq ${mapq}"
    sbatch --time=$jtime --mem-per-cpu=$jme --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
done
