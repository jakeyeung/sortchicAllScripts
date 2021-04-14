#!/bin/sh
# Jake Yeung
# 6-make_count_and_RZ_tables.sh
#  
# 2019-12-21

# WRAP UP


# WRAP UP
while [[ `qstat | grep merge | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# inmain=$wd
jsuffix="imputevarfilt.lessstringent.mapq_40"
inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jsuffix}.remerged_by_marks"
cd $inmain

jmem='32G'
jtime='6:00:00'

mapq=40
binsize=10000  # try other binsizes?
stepsize=${binsize}
# ps="/home/hub_oudenaarden/jyeung/projects/SingleCellMultiOmics.ForDev/singlecellmultiomics/bamProcessing/bamToCountTable.test.py"
ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-06-03/singlecellmultiomics/bamProcessing/bamToCountTable.py"
bl="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.lessstringent.2020-04-14.corrfilt/correlated_bins.chromofixed.bed"

outdir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_${binsize}.${jsuffix}"
[[ ! -d $outdir ]] && mkdir $outdir

for inbami in `ls -d $inmain/*.bam.bai`; do
    inbam=${inbami%.*}
    [[ ! -e $inbam ]] && echo "$inbam not found, continuing" && continue 
    bname=$(basename $inbam)
    bname=${bname%.*}

    BNAME=$outdir/${bname}.counttables.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    BNAME2=$outdir/${bname}.lhcounts.qsub
    DBASE2=$(dirname "${BNAME2}")
    [[ ! -d $DBASE2 ]] && echo "$DBASE2 not found, exiting" && exit 1

    outf=$outdir/${bname}.countTable.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam -sliding $stepsize --filterXA -minMQ $mapq -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.countTable

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps -bin $binsize -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name --dedup -blacklist $bl --r1only"
    sbatch --time=${jtime} --mem-per-cpu=${jmem} --output=${BNAME}_%j.log --ntasks=1 --cpus-per-task=1 --nodes=1 --ntasks-per-node=1 --ntasks-per-socket=1 --job-name=$prefix.countTablesTSS --wrap "$cmd"

    # outf2=$outdir/${bname}.LHcounts.csv
    # [[ -e $outf2 ]] && echo "$outf2 found, continuing" && continue
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inbam -sampleTags SM -featureTags lh -o $outf2 --dedup --filterXA -minMQ $mapq" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME2}.out -e ${BNAME2}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname.RZcounts

done

