#!/bin/sh
# Jake Yeung
# 4-run.filter_mat_by_good_cells.sh
#  
# 2019-11-08

jmem='4G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/Zebrafish/make_count_tables_TSS/filter_mat_by_good_cells.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
ouds="oud3985"

goodcells="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_ZFbonemarrow/good_cells/ZF-H3K27me3.2019-11-13.csv"

for oud in $ouds; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/${oud}/countTables_geneTSS"
    [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
    for inf in `ls -d $inmain/PZ-ChIC-ZFWKM-H3K27me3*.rds`; do
        bname=$(basename $inf)
        bname=${bname%.*}
        outname="$bname.goodcellsfilt.RData"
        BNAME=$inmain/$bname.goodcellsfilt.qsub
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        outf=$inmain/$outname
        [[ -e $outname ]] && echo "$outname found, continuing" && continue
        echo "Rscript $rs $inf $goodcells $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname
    done
done



