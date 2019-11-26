#!/bin/sh
# Jake Yeung
# 6-make_RZ_counts.sh
#  
# 2019-09-04

# sleep 3600

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all/raw_demultiplexed"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all/RZcounts"
[[ ! -d $outdir ]] && mkdir $outdir

jmem='8G'
jtime='1:00:00'

for indir in `ls -d $inmain/PZ-ChIC-ZFWKM*`; do
    bname=$(basename $indir)
    inbam=$indir/tagged/$bname.bwaMapped.tagged.bam
    [[ ! -e $inbam ]] && echo "$inbam not found, skipping" && continue
    outf=$outdir/${bname}.LH_counts.demuxbugfixed.csv

    BNAME=$outdir/${bname}.LHcounts.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py $inbam -sampleTags SM -featureTags RZ -o $outf --dedup --filterXA -minMQ 40" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname.RZcounts
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamToCountTable.py $inbam -sampleTags SM -featureTags lh -o $outf --dedup --filterXA -minMQ 40" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname.RZcounts
done

