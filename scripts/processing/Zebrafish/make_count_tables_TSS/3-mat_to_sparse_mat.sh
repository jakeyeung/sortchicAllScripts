#!/bin/sh
# Jake Yeung
# 5-mat_to_sparse_mat.sh
#  
# 2019-06-17

jmem='4G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/1-qc_filter_bins_cells_bin_matrix/lib/mat_to_sparse_mat_new_slidewin_format.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# oud="oud3910"
ouds="oud3909 oud3910"

for oud in $ouds; do
  inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/${oud}_HD_0/countTables_geneTSS"
  outmain="$inmain"
  [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
  
  # . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6
  # Rscript="/hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/Rscript"
  
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
done

