#!/bin/sh
# Jake Yeung
# 11-run_LDA_on_rds.sh
#  
# 2019-11-13

jmem='42G'
jtime='150:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

# bsize="50000_50000"
bsize="20000_20000"
# bsize="20000_10000"

prefix="B6BM_All_allmarks.2020-01-12.bsizestepsize_${bsize}"

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_B6.2019-12-16"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_redo_2019-12-13.tagged_bams_mergedbymarks/countTables_otherWinSize_NoSliding/qc_filt/qc_${bsize}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_${prefix}"  
[[ ! -d $outmain ]] && mkdir $outmain
[[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1  

for inf in `ls -d $indir/count_mat*.rds`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}"
    [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N bsize_${bsize}_$bname
    # echo "cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname"
done
