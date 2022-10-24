#!/bin/sh
# Jake Yeung
# 4-run_LDA_doubles.sh
#  
# 2021-11-21

jmem='16G'
jtime='72:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scChIX/utils/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
# binarize="FALSE"

cell="K562"
outmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/LDA_outputs/${cell}"
[[ ! -d $outmain0 ]] && mkdir $outmain0

jmarks="k4me1 k4me3 k27me3 k9me3"
for jmark in $jmarks; do 
  jname="${cell}_${jmark}"

  inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/filtered_count_tables_for_LDA/${cell}/${jname}"
  [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
  
  for inf in `ls -d ${inmain}/*.rds`; do
      prefix=$(basename $inf)
      prefix=${prefix%.*}
      [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
      echo $inf
      outmain="${outmain0}/ldaAnalysis_fripfilt_${jname}"
      [[ ! -d $outmain ]] && mkdir $outmain
      [[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1
      [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
  
      bname="${prefix}"
  
      outdir="${outmain}/lda_outputs.${bname}"
      [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
      [[ ! -d $outdir ]] && mkdir -p $outdir
  
      BNAME=$outdir/$bname
      DBASE=$(dirname "${BNAME}")
      [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
      cmd="cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname"
      sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=LDA_${bname} --wrap "$cmd"
  done

done