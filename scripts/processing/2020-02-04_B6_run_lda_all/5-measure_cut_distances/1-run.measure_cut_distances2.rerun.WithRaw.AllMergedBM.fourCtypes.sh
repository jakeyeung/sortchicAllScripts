#!/bin/sh
# Jake Yeung
# 1-run.measure_cut_distances.sh
#  
# 2020-07-09

jmem='32G'
jtime='4:00:00'
ncores=1

ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-07-17.NewCountFilters/singlecellmultiomics/bamProcessing/bamAnalyzeCutDistances.JY.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.withBlood.bams_remerged.MAPQ_40"
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster.MAPQ_40"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/DedupR1onlyNoAltHits"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/cut_distances_all_clusters2.rerun.withraw.BM.fourCtypes"
[[ ! -d $outmain ]] && mkdir $outmain

for inbam in `ls -d $indir/*.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    outdir=${outmain}/${bname}
    [[ ! -d $outdir ]] && mkdir $outdir
    BNAME=${outdir}/${bname}.sbatchlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $inbam -o $outdir"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=$ncores --nodes=1 --job-name=${bname} --wrap "$cmd"
done

