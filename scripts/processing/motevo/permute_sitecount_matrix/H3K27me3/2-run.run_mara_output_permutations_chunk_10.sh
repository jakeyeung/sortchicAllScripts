#!/bin/sh
# Jake Yeung
# 2-run.run_mara_output_permutations_chunk_10.sh
# Run 1000 chunks of 10 
# 2019-04-01

jmark="H3K27me3"

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/run_mara_output_permutations_chunk_10.sh"

[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

chunkstart=2
chunkend=500

jmem='16G'
jtime='1:00:00'
by="row"

nohupdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${jmark}/nohups"
[[ ! -d $nohupdir ]] && mkdir $nohupdir


for chunk in $(seq $chunkstart $chunkend); do
# for chunk in $(seq 1); do
    BNAME=${nohupdir}/nohup_${chunk}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo "bash $bs $chunk $by" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 10 -m beas -M j.yeung@hubrecht.eu
    # echo "bash $bs $chunk $by $jmark" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 10
    echo "bash $bs $chunk $by $jmark" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 10 
    # bash $bs $chunk $by $jmark
    # exit 0
done
