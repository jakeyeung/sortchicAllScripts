#!/bin/sh
# Jake Yeung
# run.split_bam_by_cluster.sh
#  
# 2020-01-09

jmem='32G'
jtime='12:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/split_bam_by_cluster.py"
jmarks="H3K4me3 H3K27me3 H3K9me3"
jwin="50000_25000"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks"
annotdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_clustering_BM"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $annotdir ]] && echo "$annotdir not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_redo_2019-12-13.cluster_analysis/bams_split_by_cluster"

[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1  # make dir before hand

for jmark in $jmarks; do
    for inf in `ls -d ${indir}/${jmark}*.bam`; do
        bname=$(basename $inf)
        bname=${bname%.*}
        annotfile=${annotdir}/clusters_BM_All_merged_${jmark}_${jwin}.txt
        [[ ! -e $annotfile ]] && echo "$annotfile not found, exiting" && exit 1
        BNAME=$outdir/$bname
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq 40 --add_chr_prefix" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N BM_${jmark}
    done
done

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/H3K4me3-BM_Linneg_SC-merged.tagged.bam"
# annotfile="/hpc/hub_oudenaarden/jyeung/data/intestinal_scchic/from_rstudioiserver/pdfs_clustering/clusters_BM_All_merged_H3K4me3_20000_10000.txt"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_redo_2019-12-13.split_bams_by_clusters"

# . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -infile $inf -annotfile $annotfile -outdir $outdir -mapq 40 --add_chr_prefix

