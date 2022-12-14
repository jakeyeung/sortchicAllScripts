#!/bin/sh
# Jake Yeung
# 8-run.run_GREAT_downstream.sh
# Run run_GREAT_downstream.R on cluster
# 2019-01-06

ncores=5

jmem='32G'
jtime='4:00:00'

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/run_GREAT_downstream_HS.R"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks_0.5_1000/lda_outputs.meanfilt_1.cellmin_1000.cellmax_50000"
# maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks_0.5_1000/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000"
# maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks_0.5_1000/lda_outputs.meanfilt_0.32.cellmin_1000.cellmax_50000"
jdist=1000
jmeanfilt=1
jcellmin=100
jcellmax=500000
meth="hiddenDomains"
binarize="TRUE"
cell="K562"

# maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysis${meth}_${jdist}/lda_outputs.meanfilt_${jmeanfilt}.cellmin_${jcellmin}.cellmax_${jcellmax}.binarize.${binarize}"

maindir="/hpc/hub_oudenaarden/jyeung/data/scChiC/LDA_out_${cell}/ldaAnalysishiddenDomains_${jdist}_round2/lda_outputs.hiddenDomains.meanfilt_${jmeanfilt}.cellmin_${jcellmin}.cellmax_${jcellmax}.binarize.${binarize}"
[[ ! -d $maindir ]] && echo "$maindir not found, exiting" && exit 1

chainpath="/hpc/hub_oudenaarden/jyeung/data/databases/chainfiles/hg38ToHg19.over.chain"
[[ ! -e $chainpath ]] && echo "$chainpath not found, exiting" && exit 1

outdir=$maindir/downstream
[[ ! -d $outdir ]] && mkdir $outdir

for inf in $(ls -d $maindir/*.Robj); do
    bname=$(basename $inf)
    bname2=${bname%.*}  # remove .Robj extension
    outf=$outdir/$bname2.GREAT.Robj
    [[ -e $outf ]] && echo "$outf found, skipping" && continue
    BNAME=$outdir/$bname2.GREATjob
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo "cd $wd; Rscript $rs $inf $outf $ncores $chainpath" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu
    # echo "cd $wd; Rscript $rs $inf $outf $ncores"
    # cd $wd; Rscript $rs $inf $outf $ncores
    # exit 0
done
