#!/bin/sh
# Jake Yeung
# 10-run.filter_bins_by_chromosomes.sh
#  
# 2020-04-08

jmem='16G'
jtime='2:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-04-03_ZebrafishAll/processing/filter_bins_by_chromosomes.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000"
# inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000.filtered_by_counts_TAfrac_var"
inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.winsize_50000.chromofilt"
[[ ! -d $outmain ]] && mkdir $outmain

chromos="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chr23 chr24 chr25"

# for inf in `ls -d $inmain/*.rds`; do 
for inf in `ls -d $inmain/*.countTable.csv`; do 
  bname=$(basename $inf)
  bname=${bname%.*}
  bnameout=${bname}.chrfilt.rds
  
  outf=${outmain}/${bnameout}

  BNAME=$outmain/qsub_${bname}
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

  echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf -chromoskeep $chromos" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $bname
done
