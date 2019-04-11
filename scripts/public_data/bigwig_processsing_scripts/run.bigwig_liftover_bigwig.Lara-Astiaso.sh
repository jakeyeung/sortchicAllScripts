#!/bin/sh
# Jake Yeung
# run.bigwig_liftover_bigwig.Lara-Astiaso.sh
# bigwig liftover bigwig  
# 2019-03-28

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/public_data/bigwig_liftover_bigwig.ForServer.sh"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_GenomeRes_2014/bigwig_mm9"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/bigwig_mm10"
[[ ! -d $outdir ]] && mkdir $outdir

n=0
maxjobs=12
for w in `ls -d $indir/*.bigWig`; do
                # echo "bash $rs $w $outdir&"
                bash $rs $w $outdir&
                if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
                        # define maxjobs and n using maxjobsn skeleton
                    wait # wait until all have finished (not optimal, but most times good enough)
                    echo $n wait
                fi
done
wait
