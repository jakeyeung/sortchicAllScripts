#!/bin/sh
# Jake Yeung
# 7-make_count_tables.sh
#  
# 2019-09-04


inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks"
outdir=$inmain/countTables_otherWinSize
[[ ! -d $outdir ]] && mkdir $outdir

jmem='64G'
jtime='4:00:00'

mapq=40
# stepsize=20000
# binsize=10000
bsizes="5000 10000 20000 50000"

for bsize in $bsizes; do
    stepsize=$(expr $bsize / 2)
    # stepsize=$(echo $((x / y)))
    # inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_mergedAll.retag/H3K4me3-BM_SC-merged.tagged.retagged.bam"
    for inbam in `ls -d $inmain/*.bam`; do
        [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
        bname=$(basename $inbam)
        bname=${bname%.*}
        outf=$outdir/${bname}.bsize_${bsize}.step_${stepsize}.countTable.demuxbugfixed.csv
        [[ -e $outf ]] && echo "$outf found, continuing" && continue

        BNAME=$outdir/${bname}.counttables.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $bsize -binTag DS --dedup" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.$bsize.countTable
        . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $bsize -binTag DS --dedup
        exit 0
        # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $bsize -binTag DS --dedup&
    done
done
wait

