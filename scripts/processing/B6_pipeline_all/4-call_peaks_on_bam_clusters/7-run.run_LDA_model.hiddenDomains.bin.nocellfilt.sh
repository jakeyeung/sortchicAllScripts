#!/bin/sh
# Jake Yeung
# 2b-run.run_LDA_model.H3K27me3.BroadPeaks.sh
#  
# 2018-12-29

n=0
maxjobs=1

# if not binarize
jmem='60G'
jtime='32:00:00'
# if binarize
# jmem='48G'
# jtime='60:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/lib/run_LDA_model.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
# marks="H3K4me3"
marks="H3K27me3"
# marks="H3K4me3"
# marks="H3K27me3"
# marks="H3K9me3"
mindist="1000"
cell="BM"

K=20  # kind of useless parameter
ncores=1
# topics="50"
topics="50"  # maybe Tal1 shows up higher?
# topics="30"  # maybe Tal1 shows up lower?
topicsName=`echo $topics | sed 's/,/_/g'`
tunemodels="TRUE"
binarize="TRUE"
cellmin="0"
cellmax="9999999"
meanmax="10"  # meanmax = 1 keeps the Erdr1 and Mid1 genes, which are probably real? But there may be some weird peaks as well that are skewing results.

suffix="build95_B6.cells_from_bin_analysis"
# suffix2="GeneTSS.Dedup"
# suffix2="GeneTSS.Dedup.RbindHiddenDomains"
suffix2="CorrPeakFilt"
# tssdist=20000

# marks="H3K27me3"
for mark in $marks; do
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.${mindist}_${suffix}/PZ-${cell}-${mark}.merged.NoCountThres.hiddenDomains.Robj"
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.1000_build95.withchr.cells_from_bin_analysis/PZ-BM-${mark}.merged.NoCountThres.GeneTSS.Dedup.RbindHiddenDomains.Robj"
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.1000_build95.withchr_${tssdist}.cells_from_bin_analysis/PZ-BM-${mark}.merged.NoCountThres.${suffix2}.Robj"
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.${mindist}_${suffix}/CorrPeakFilt/B6-${mark}.merged.NoCountThres.hiddenDomains.${suffix2}.Robj"
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.${mindist}_${suffix}/CorrPeakFilt/PZ-BM-${mark}.merged.NoCountThres.hiddenDomains.${suffix2}.Robj"

    # relax assumptions to capture more H3K4me3 and H3K9me3 cells?
    # outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysis_${suffix2}_${mindist}_${suffix}_${tssdist}_2019-04-20/lda_outputs.meanfilt_${meanmax}.cellmin_${cellmin}.cellmax_${cellmax}.binarize.${binarize}"
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysis_${suffix2}_${mindist}_${suffix}/lda_outputs.meanfilt_${meanmax}.cellmin_${cellmin}.cellmax_${cellmax}.binarize.${binarize}"
    [[ ! -d $outdir ]] && mkdir -p $outdir

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

    bname=$(basename $inf)
    bname=${bname%%.*}.CountThres0.K-${topicsName}
    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo $bname
    # echo "cd $workdir; Rscript $rs $inf $outdir $K $topics $tunemodels $meanmax $cellmin $cellmax $binarize $bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu
    cd $workdir; Rscript $rs $inf $outdir $K $topics $tunemodels $meanmax $cellmin $cellmax $binarize $bname&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
