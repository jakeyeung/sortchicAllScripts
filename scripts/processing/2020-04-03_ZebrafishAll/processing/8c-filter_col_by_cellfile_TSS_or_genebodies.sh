#!/bin/sh
# Jake Yeung
# 8c-filter_row_by_cellfile_TSS_or_genebodies.sh
#  
# 2020-04-14

# WRAP UP
while [[ `qstat | grep filt_ | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

jmem='16G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/filter_good_cells_by_cellfile.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.TSS.winsize_50000"
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.genebody"
outdir=${indir}
celldir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.lessstringent.2020-04-14"

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for mark in $marks; do
    # inf=${indir}/"PZ-ChIC-ZF_${mark}_2020-04-07.countTable.TSS.rds"
    inf=${indir}/"PZ-ChIC-ZF_${mark}_2020-04-07.countTable.genebody.rds"
    fname=$(basename $inf)
    fname=${fname%.*}
    cellfile=${celldir}/"count_mat.${mark}.countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.colnames"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ ! -e $cellfile ]] && echo "$cellfile not found, exiting" && exit 1
    outf="${outdir}/${fname}.cellfilt.rds"

    BNAME=$outdir/qsubCellfilt_${fname}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $cellfile $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N filt_col_$mark 
done

