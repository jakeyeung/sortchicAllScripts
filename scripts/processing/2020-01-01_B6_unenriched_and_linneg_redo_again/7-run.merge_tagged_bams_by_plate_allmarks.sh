#!/bin/sh
# Jake Yeung
# 7-merge_tagged_bams_by_plate_H3K4me1.sh
# 
# 2019-12-14

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_unenriched_and_linneg_redo_again/7-merge_tagged_bams_by_plate_allmarks.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

# WRAP UP
while [[ `qstat | wc -l` > 1 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='64G'
jtime='12:00:00'
BNAME="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13/tmpdir/mergebams"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo "bash $bs" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N mergebams -m beas -M j.yeung@hubrecht.eu

