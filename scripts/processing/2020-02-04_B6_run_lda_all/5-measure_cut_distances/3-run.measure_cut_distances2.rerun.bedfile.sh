#!/bin/sh
# Jake Yeung
# 1-run.measure_cut_distances.sh
#  
# 2020-07-09

jmem='16G'
jtime='2:00:00'
ncores=1

# ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-06-30/singlecellmultiomics/bamProcessing/bamAnalyzeCutDistances.py"
ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-07-17.NewCountFilters/singlecellmultiomics/bamProcessing/bamAnalyzeCutDistances.py"

# inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/H3K9me3-BM_AllMerged.Neutrophils_topic9.sorted.bam"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

# ctypes="Eryth Bcell Granu HSPCs"
ctypes="Erythroblast HSCs Bcell Neutrophil"
# ctype="Bcell"
bedmain="/hpc/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/bedannotations/MouseBMFromRNAseq"
dist="50000"

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/cut_distances_all_clusters2.rerun.bedfile_winsize_${dist}"
[[ ! -d $outmain ]] && mkdir $outmain

for ctype in $ctypes; do
    bedfile="${bedmain}/MouseBM_TSS.${ctype}.bsize_2.bed"
    [[ ! -e $bedfile ]] && echo "$bedfile not found, exiting" && exit 1

    for inbam in `ls -d $indir/*.bam`; do
        bname=$(basename $inbam)
        bname=${bname%.*}
        outdir=${outmain}/geneset_${ctype}_${bname}
        [[ ! -d $outdir ]] && mkdir $outdir
        BNAME=${outdir}/${bname}.sbatchlog
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $inbam -o $outdir -regions $bedfile -region_radius ${dist}"
        # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $inbam -o $outdir -regions $bedfile
        # sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=8 --nodes=1 --ntasks-per-node=8 --ntasks-per-socket=8 --job-name=${bname} --wrap "$cmd"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=$ncores --nodes=1 --job-name=${bname} --wrap "$cmd"
        # exit 0
    done
done
