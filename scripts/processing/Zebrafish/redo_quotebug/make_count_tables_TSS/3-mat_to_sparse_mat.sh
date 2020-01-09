#!/bin/sh
# Jake Yeung
# 5-mat_to_sparse_mat.sh
#  
# 2019-06-17

jmem='4G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/1-qc_filter_bins_cells_bin_matrix/lib/mat_to_sparse_mat_new_slidewin_format.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag/countTables_TSS"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outmain=$inmain/rds_outputs
[[ ! -d $outmain ]] && mkdir $outmain

for inf in `ls -d $inmain/*.csv.gz`; do
    echo $inf
    bname=$(basename $inf)
    bname=${bname%.*}  # remove .gz
    bname=${bname%.*}  # remove .csv
    BNAME=$outmain/qsub_sparsemat_$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    outf="$outmain/${bname}.rds"
    # [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf&
done
wait
