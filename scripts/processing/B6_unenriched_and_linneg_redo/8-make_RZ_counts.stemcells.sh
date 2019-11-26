#!/bin/sh
# Jake Yeung
# 8-make_RZ_counts.linneg.sh
#  
# 2019-11-21

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6StemCells/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6StemCells/RZcounts"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='8G'
jtime='2:00:00'

for indir in `ls -d $inmain/*BM*`; do
    bname=$(basename $indir)
    inbam=$indir/tagged/${bname}.bwaMapped.sorted.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, skipping" && continue
    outf=$outdir/${bname}.LH_counts.demuxbugfixed.csv
    [[ ! -e $outf ]] && echo "$outf not found, continuing" && continue

    BNAME=$outdir/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py $inbam -sampleTags SM -featureTags lh -o $outf --dedup --filterXA -minMQ 40" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname.LHcounts
done

