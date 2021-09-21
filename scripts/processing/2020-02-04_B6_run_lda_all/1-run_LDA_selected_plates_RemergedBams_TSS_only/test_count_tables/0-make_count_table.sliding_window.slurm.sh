#!/bin/sh
# Jake Yeung
# 1-bam_to_count_table_TSS.sh
#  
# 2019-06-23

jmem='32G'
jtime='02:00:00'

mapq=40
# dist=10000
prefixs="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

outsuffix=""

ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.WithBlacklist.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_remerged_by_cluster.MAPQ_40"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1


bl="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.bed"
[[ ! -e $bl ]] && echo "$bl not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow${outsuffix}.slurm.noR2.AddOne"
[[ ! -d $outdir ]] && mkdir $outdir

# bsize=10000
# stepsize=10000
size=10000
bsize=$size
stepsize=$size
dist=$size

for prefix in $prefixs; do
    inbam="${inmain}/${prefix}*.bam"
    outf=${outdir}/${prefix}.mapq_${mapq}.SlidingWindow_${dist}.blfiltered.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$outdir/${prefix}.${dist}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # python $ps -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $bsize -binTag DS --dedup -blacklist $blfile

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -bin $bsize -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -blacklist $bl" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu  -N $prefix.countTablesTSS
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -bin $bsize -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -blacklist $bl" | sbatch --time=${jtime} --mem-per-cpu=${jmem} --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=$prefix.countTablesTSS
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -bin $bsize -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -blacklist $bl"
    sbatch --time=${jtime} --mem-per-cpu=${jmem} --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=$prefix.countTablesTSS --wrap "$cmd"
done
wait
