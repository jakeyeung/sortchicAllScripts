#!/bin/sh
# Jake Yeung
# 7-make_count_tables.sh
#  
# 2019-09-04

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_linneg_K4me1.tagged_bams_mergedbymarks"
outdir=$inmain/countTables
[[ ! -d $outdir ]] && mkdir $outdir

jmem='64G'
jtime='6:00:00'

mapq=40
stepsize=50000
binsize=50000

for inbam in `ls -d $inmain/*.tagged.bam`; do
    bname=$(basename $inbam)
    bname=${bname%.*}
    outf=$outdir/${bname}.countTable.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$outdir/${bname}.counttables.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname.countTable
done

