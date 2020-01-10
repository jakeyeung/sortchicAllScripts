#!/bin/sh
# Jake Yeung
# 1-run.mat_to_mm.sh
#  
# 2019-12-16

jmem='4G'
jtime='1:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/sparse_mat_to_mm.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all.retag/countTables_TSS/rds_outputs/good_cells_filt"
outmain="$inmain/for_gensim"
[[ ! -d $outmain ]] && mkdir $outmain

for inf in `ls -d $inmain/*.RData`; do
    bbase=$(basename $inf)
    bbase=${bbase%.*}
    outf=$outmain/$bbase.mm
    BNAME=$outmain/$bbase.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bbase
done
