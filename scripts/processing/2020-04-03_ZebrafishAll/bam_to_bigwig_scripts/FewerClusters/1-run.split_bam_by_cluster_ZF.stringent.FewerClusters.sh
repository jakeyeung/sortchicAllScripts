#!/bin/sh
# Jake Yeung
# run.split_bam_by_cluster.sh
#  
# 2020-01-09

jmem='16G'
jtime='12:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# jmarks="H3K4me3"
mapq=40

jdate="2020-04-23"
jprefix="imputevarfilt.lessstringent"
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks"
# annotdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain"
# annotdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.${jdate}.${jprefix}"
annotdir="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/cell_to_cluster_tables_merged"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $annotdir ]] && echo "$annotdir not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jprefix}.mapq_${mapq}.FewerClusters.2020-06-17"
[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1  # make dir before hand

for jmark in $jmarks; do
    bname="PZ-ChIC-ZF_${jmark}_2020-04-07"
    inf=${indir}/${bname}.bam
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    # annotname="cell_to_cluster_table.${jmark}.2020-04-13.txt"
    # annotname="ZF_LDA_output.${jmark}.keepn_150.final.ClusterTables.txt"
    annotname="WKM_cell_to_clusters.${jmark}.txt"
    annotfile=${annotdir}/${annotname}
    [[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq ${mapq} --add_chr_prefix"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${bname}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=${bname} --wrap "$cmd"
done

