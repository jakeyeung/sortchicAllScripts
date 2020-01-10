#!/bin/sh
# Jake Yeung
# 7b-check_bams.sh
#  
# 2019-12-15

# check bams

jmem='1G'
jtime='2:00:00'

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks"

for bfull in `ls -d $inf/*.bam.bai`; do
    # echo $bfull
    bamin=${bfull%.*}
    bmain=${bamin%.*}
    bmainbase=$(basename $bmain)
    bout=$bmain.check.out
    # echo $bamin
    BNAME=$bmain.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; bamTabulator.py $bamin SM > $bout" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N check.$bmainbase
done
