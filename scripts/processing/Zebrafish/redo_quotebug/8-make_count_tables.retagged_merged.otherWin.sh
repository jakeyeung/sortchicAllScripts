#!/bin/sh
# Jake Yeung
# 8-make_count_tables.retagged_merged.sh 
# Retag to get RC tag, merge across plates, then make tables
#  
# 2019-09-04

# sleep 3600

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag/merge_by_mark"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag/countTables_otherWinSize"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='12G'
jtime='6:00:00'

mapq=40
# stepsize=20000
# binsize=100000
bsizes="10000 20000 50000"

inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag/merge_by_mark/H3K4me1-WKM_CD41-merged.tagged.retagged.bam"
# for inbam in `ls -d $inmain/*.retagged.bam`; do  # zebrafish only 
for binsize in $bsizes; do
    stepsize=$(expr $binsize / 2)
    bname=$(basename $inbam)
    bname=${bname%.*}
    outf=$outdir/${bname}.bsize_${binsize}.stepsize_${stepsize}.countTable.demuxbugfixed_mergeplates.csv
    [[ -e $outf ]] && echo "$outf found, continuing" && continue

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py -sliding $stepsize --filterXA -minMQ $mapq $inbam -o $outf -sampleTags SM -joinedFeatureTags reference_name -bin $binsize -binTag DS --dedup&
done
wait

