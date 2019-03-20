#!/bin/sh
# Jake Yeung
# 8-run.run_GREAT_downstream.sh
# Run run_GREAT_downstream.R on cluster
# 2019-01-06

ncores=10

jmem='32G'
jtime='4:00:00'

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/run_GREAT_downstream.R"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1


jbin="TRUE"
# maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.FALSE.no_filt/lda_out_meanfilt.BM-H3K27me3.CountThres0.K-15_20_25_30.Robj"
maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.${jbin}.no_filt"

[[ ! -d $maindir ]] && echo "$maindir not found, exiting" && exit 1

outdir=$maindir/downstream
thres=0.96
[[ ! -d $outdir ]] && mkdir $outdir

for inf in $(ls -d $maindir/*.Robj); do
    bname=$(basename $inf)
    bname2=${bname%.*}  # remove .Robj extension
    outf=$outdir/$bname2.GREAT.${thres}.Robj
    [[ -e $outf ]] && echo "$outf found, skipping" && continue
    BNAME=$outdir/$bname2.GREATjob.${thres}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo "cd $wd; Rscript $rs $inf $outf $ncores $thres" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu
done
