#!/bin/sh
# Jake Yeung
# 6-add_genename_to_rows.sh
# Gene names in rows were lost during mat to sparse mat 
# add it back!
# 2019-11-11

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/Zebrafish/make_count_tables_TSS/add_genename_to_rows.R"
[[ ! -e $rs  ]] && echo "$rs  not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_oud3910_merged/countTables_geneTSS/merged"
annotmain="/hpc/hub_oudenaarden/jyeung/data/databases/gene_tss"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_oud3910_merged/countTables_geneTSS/mergedRownamesReadded"
[[ ! -d $outmain ]] && mkdir $outmain

for inf in `ls -d $inmain/*.RData`; do
  infbase=$(basename $inf)
  winsize=`echo "$infbase" | cut -d"." -f2 | cut -d"_" -f2`
  echo $infbase
  annot=${annotmain}/"gene_tss.winsize_${winsize}.species_drerio.nochr.bed"
  [[ ! -e $annot ]] && echo "$annot not found, exiting" && exit 1

  outf="${outmain}/${infbase}"
  # echo "Rscript $rs $inf $annot $outf"
  Rscript $rs $inf $annot $outf
done


