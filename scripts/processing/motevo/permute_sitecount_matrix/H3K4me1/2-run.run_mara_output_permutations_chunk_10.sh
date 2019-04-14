#!/bin/sh
# Jake Yeung
# 2-run.run_mara_output_permutations_chunk_10.sh
# Run 1000 chunks of 10 
# 2019-04-01

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/run_mara_output_permutations_chunk_10.sh"

chunkstart=1
chunkend=500

jmem='24G'
jtime='1:00:00'

nohupdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_H3K4me1/nohups"

by="row"

for chunk in $(seq $chunkstart $chunkend); do
# for chunk in $(seq 1); do
    BNAME=${nohupdir}/nohup_${chunk}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo "bash $bs $chunk $by" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 10 -m beas -M j.yeung@hubrecht.eu
    echo "bash $bs $chunk $by" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 10
done
