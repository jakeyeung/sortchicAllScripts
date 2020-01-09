#!/bin/sh
# Jake Yeung
# 4-run.filter_mat_by_good_cells.sh
#  
# 2019-11-08

jmem='4G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/Zebrafish/make_count_tables_TSS/filter_mat_by_good_cells.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
goodcells="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_ZF/good_cells/good_cells_all_merged.csv"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag/countTables_TSS/rds_outputs"
outdir=$inmain/good_cells_filt
[[ ! -d $outdir ]] && mkdir $outdir

[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
for inf in `ls -d $inmain/*.rds`; do
    bname=$(basename $inf)
    bname=${bname%.*}
    outname="$bname.goodcellsfilt.RData"
    BNAME=$outdir/$bname.goodcellsfilt.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    outf=$outdir/$outname
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $goodcells $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname
    . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $goodcells $outf&
done
wait




