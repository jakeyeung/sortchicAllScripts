#!/bin/sh
# Jake Yeung
# 8-run.run_metacell.sh
# Run metacell  
# 2018-12-20

jmem='128G'
jtime='8:00:00'

rc="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/run_metacell.R"

[[ ! -e $rc ]] && echo "$rc not found, exiting" && exit 1

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats/PZ-BM-H3K4me1.merged.NoCountThres.subset.mat"  # subset to iron out bugs
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats/PZ-BM-H3K4me1.merged.NoCountThres.mat"  # prime time

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/mc_outputs_full"
[[ -d $outdir ]] && echo "$outdir must not exist" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

# jvar=0.05

BNAME="$outdir/PZ-BM-H3K4me1.merged.NoCountThres"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

# echo "Rscript $rc $inf $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu 
# echo "Rscript $rc $inf $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu 
echo "Rscript $rc $inf $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu 
