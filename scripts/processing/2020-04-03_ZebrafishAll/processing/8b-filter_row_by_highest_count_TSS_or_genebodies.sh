#!/bin/sh
# Jake Yeung
# 8b-filter_row_by_highest_count_TSS_or_genebodies.sh
#  
# 2020-04-14

jmem='16G'
jtime='1:00:00'


rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/processing/filter_row_by_highest_count.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.TSS.winsize_50000"
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.genebody"
outdir="${indir}"

for inf in `ls -d $indir/*.genebody.csv`; do
    fname=$(basename $inf)
    fname=${fname%.*}
    outf=${outdir}/${fname}.rds
    BNAME=$outdir/qsubBinfilt_${fname}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf -format bed" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N filt_$fname
done
