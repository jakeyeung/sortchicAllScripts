#!/bin/sh
# Jake Yeung
# 5-mat_to_sparse_mat.sh
#  
# 2019-06-17

jmem='4G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/mat_to_sparse_mat_tss_format.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1


inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataZF_all/countTables_geneTSS"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outmain="$inmain"


for inf in `ls -d $inmain/PZ-ChIC*.csv`; do
  echo $inf
  bname=$(basename $inf)
  bname=${bname%.*}
  BNAME=$outmain/qsub_sparsemat_$bname
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  outf="$outmain/${bname}.rds"
  # [[ -e $outf ]] && echo "$outf found, continuing" && continue
  echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -N $bname
done
wait

