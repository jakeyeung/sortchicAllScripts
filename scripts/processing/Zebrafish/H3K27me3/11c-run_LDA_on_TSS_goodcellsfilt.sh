#!/bin/sh
# Jake Yeung
# 11c-run_LDA_on_TSS_goodcellsfilt.sh
#  
# 2019-11-13

jmem='24G'
jtime='24:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

K=20  # kind of useless parameter
ncores=3
topics="30,35,50"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

prefix="ZFbonemarrow"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3985/countTables_geneTSS/merged"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_${prefix}"  # add to existing directory
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

for inf in `ls -d $indir/PZ-ChIC-ZFWKM-H3K27me3*merged.RData`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}"
    # [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N $bname
    # echo "cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname"
done
