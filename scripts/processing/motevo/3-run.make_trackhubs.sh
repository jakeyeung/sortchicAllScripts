#!/bin/sh
# Jake Yeung
# 3-run.make_trackhubs.sh
#  
# 2019-02-01

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py2

script="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/make_trackhubs.py"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/bigbeds"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/trackhub_files/H3K4me1_peaks"

[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

python $script $indir $outdir --render --suffix "H3K4me1"
