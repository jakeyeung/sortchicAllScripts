#!/bin/sh
# Jake Yeung
# 2-run.run_mara_output_permutations_chunk_10.sh
# Run 1000 chunks of 10 
# 2019-04-01

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/run_mara_output_permutations_chunk_10_general.sh"

chunkstart=1
chunkend=200

jmem='108G'
jtime='3:00:00'

nohupdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_H3K4me1/nohups_CorrPeakFilt"
[[ ! -d $nohupdir ]] && mkdir $nohupdir

by="row"

mark="H3K4me1"  # USELESS PARAMETER
exprsrds="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_H3K4me1/mara_input/exprs_mat.rds"
sitecountsrds="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_H3K4me1/mara_input/sitecounts_mat.rds"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_H3K4me1/mara_output"

[[ ! -e $exprsrds ]] && echo "$exprsrds not found, exiting" && exit 1
[[ ! -e $sitecountsrds ]] && echo "$sitecountsrds not found, exiting" && exit 1

for chunk in $(seq $chunkstart $chunkend); do
# for chunk in $(seq 1); do
    echo $chunk
    BNAME=${nohupdir}/nohup_${chunk}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo "bash $bs $chunk $by $mark $exprsrds $sitecountsrds $outmain" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 10 -m beas -M j.yeung@hubrecht.eu
    # echo "bash $bs $chunk $by $mark $exprsrds $sitecountsrds $outmain" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
    # exit 0
done
