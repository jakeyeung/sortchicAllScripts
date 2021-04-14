#!/bin/sh
# Jake Yeung
# 4-convert_txt_to_rds_count_tables.sh
#  
# 2020-02-15

jmem='64G'
jtime='0:30:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/convert_txt_to_rds_count_table.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/count_tables_from_hiddenDomains.NewCountFilters"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir=${indir}.rds_format
[[ ! -d $outdir ]] && mkdir $outdir
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/count_tables_from_hiddenDomains.rds_format"

for inf in `ls -d $indir/*.txt`; do
    bname=$(basename $inf)
    bname=${bname%.*}
    BNAME=${outdir}/${bname}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    outf=${outdir}/${bname}.rds
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf" |  qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N Convert_${bname}
done
