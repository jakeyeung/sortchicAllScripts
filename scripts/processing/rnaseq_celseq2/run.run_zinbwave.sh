#!/bin/sh
# Jake Yeung
# run.run_zinbwave.sh
# Run run_zinbwave.R 
# 2019-01-07

jmem='16G'
jtime='2:00:00'

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/run_zinbwave.R"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_dropbox/BM-celseq2/PZ-BM-celseq2-countmatrix.dat"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/rnaseq_output/PZ-BM-celseq2-countmatrix.zinbwave.Robj"

outdir=$(dirname $outf)
BNAME=$outdir/PZ-BM-celseq2-countmatrix.zinbwave
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

ncores=16
k=20
echo "cd $wd; Rscript $rs $inf $outf $k $ncores" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores

