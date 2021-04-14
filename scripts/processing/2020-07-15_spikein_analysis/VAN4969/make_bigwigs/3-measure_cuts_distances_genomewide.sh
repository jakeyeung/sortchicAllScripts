#!/bin/sh
# Jake Yeung
# 1-run.measure_cut_distances.sh
#  
# 2020-07-09

jmem='32G'
jtime='6:00:00'
ncores=4

ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-07-17.NewCountFilters/singlecellmultiomics/bamProcessing/bamAnalyzeCutDistances.py"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/bams_split_by_clusters"

for indir in `ls -d $inmain/K562*`; do
    outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/cut_distances_all_clusters.K562_spikeins"
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
        # sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=8 --nodes=1 --ntasks-per-node=8 --ntasks-per-socket=8 --job-name=${bname} --wrap "$cmd"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=$ncores --nodes=1 --job-name=${bname} --wrap "$cmd"
    done
done

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40"
# [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
