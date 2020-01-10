#!/bin/sh
# Jake Yeung
# 3-plot_likelihood_iterations.sh
#  
# 2019-12-23

jmem='8G'
jtime='1:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/plot_likelihood_across_iterations_gensim.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_python_outputs"
outdir=${inmain}/plots
[[ ! -d $outdir ]] && mkdir $outdir



for inf in `ls -d $inmain/*.err`; do
    bname=$(basename $inf)
    bname=${bname%.*}
    BNAME=$outdir/$bname.qsub.plotlikelihood
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    outf=$outdir/${bname}_plots
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; python $ps $inf $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1
done

