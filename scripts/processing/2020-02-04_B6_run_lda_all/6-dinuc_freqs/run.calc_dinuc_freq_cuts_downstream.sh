#!/bin/sh
# Jake Yeung
# run.calc_dinuc_freq_cuts_downstream.sh
#  
# 2020-08-05

jmem='16G'
jtime='1:00:00'

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2020-02-04_B6_run_lda_all/6-dinuc_freqs/calc_dinuc_freq_cuts_downstream.py"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/again2"
outdir="${indir}/freq_cuts_downstream2"
[[ ! -d $outdir ]] && mkdir $outdir

for f in `ls -d $indir/*.csv.gz`; do
    bname=$(basename $f)
    bname=${bname%.*}.dinucfreq.csv

    BNAME=$outdir/$bname.sbatchlog
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate scmo4_NewCountFilters; python $ps $f $outdir/$bname"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=$bname --wrap "$cmd"
done
