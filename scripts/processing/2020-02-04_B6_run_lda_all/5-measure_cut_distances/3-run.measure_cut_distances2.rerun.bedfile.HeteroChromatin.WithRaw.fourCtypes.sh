#!/bin/sh
# Jake Yeung
# 1-run.measure_cut_distances.sh
#  
# 2020-07-09

jmem='16G'
jtime='2:00:00'
ncores=1

ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-07-17.NewCountFilters/singlecellmultiomics/bamProcessing/bamAnalyzeCutDistances.JY.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1


dist="50000"
# topn="2000"
# topn="5000"
# bedmain="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBMFromTopics.${topn}"
# [[ ! -d $bedmain ]] && echo "$bedmain not found, exiting" && exit 1
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/cut_distances_all_clusters2.rerun.bedfile_Heterochromatin_winsize_${dist}.heterochromatin.withraw.fourCtypes"
[[ ! -d $outmain ]] && mkdir $outmain

# ctypes="Eryth Bcell Granu HSPCs"
# for ctype in $ctypes; do
    # echo $ctype
    # bedfile="${bedmain}/MouseBM_TSS_FromTopics.${ctype}.bsize_2.bed"

    bedfile="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/offsets_hetero_and_totalcuts/MouseBM_HeteroTotalCounts_50kb_bins.2020-06-20.chrfilt.bed"
    [[ ! -e $bedfile ]] && echo "$bedfile not found, exiting" && exit 1

    for inbam in `ls -d $indir/*.bam`; do
        echo $inbam
        bname=$(basename $inbam)
        bname=${bname%.*}
        outdir="${outmain}/heterochromatin_${bname}"
        [[ ! -d $outdir ]] && mkdir $outdir
        BNAME=${outdir}/${bname}.sbatchlog
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $inbam -o $outdir -regions $bedfile -region_radius ${dist}"
        # sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=8 --nodes=1 --ntasks-per-node=8 --ntasks-per-socket=8 --job-name=${bname} --wrap "$cmd"

        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=$ncores --nodes=1 --job-name=${bname} --wrap "$cmd"
    done
# done
