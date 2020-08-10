#!/bin/sh
# Jake Yeung
# 1-run.measure_cut_distances.sh
#  
# 2020-07-09

jmem='16G'
jtime='2:00:00'
ncores=4

ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-06-30/singlecellmultiomics/bamProcessing/bamAnalyzeCutDistances.py"

# inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/H3K9me3-BM_AllMerged.Neutrophils_topic9.sorted.bam"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
# outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/cut_distances_all_clusters"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/CutDistances"
[[ ! -d $outmain ]] && mkdir $outmain

for inbam in `ls -d $indir/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    outdir=${outmain}/${bname}
    [[ ! -d $outdir ]] && mkdir $outdir
    BNAME=${outdir}/${bname}.sbatchlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo3; python $ps $inbam -o $outdir"
    # sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=8 --nodes=1 --ntasks-per-node=8 --ntasks-per-socket=8 --job-name=${bname} --wrap "$cmd"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=$ncores --nodes=1 --job-name=${bname} --wrap "$cmd"
done

