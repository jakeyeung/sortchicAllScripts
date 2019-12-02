#!/bin/sh
# Jake Yeung
# 2-project_new_samples_on_LDA.nobin.CD41plus.sh
#  
# 2019-11-04

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3"

binarize="TRUE"
Kchoose="30"

jmem='32G'
jtime='24:00:00'

for jmark in $jmarks; do

    inlda="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-ZFWKM-${jmark}_pcutoff_0.CountThres0.K-30_35_50.binarize.${binarize}/lda_out_meanfilt.ZF-ZFWKM-${jmark}_pcutoff_0.CountThres0.K-30_35_50.OutObjs.K_${Kchoose}.rds"
    inbase=$(dirname $inlda)
    inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_ZFbonemarrow/ZF-ZFWKMCD41plus-${jmark}_pcutoff_0.95_binfilt_cellfilt.2019-11-03.RData"

    [[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
    [[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

    outdir=${inbase}/projections
    [[ ! -d $outdir ]] && mkdir $outdir

    inldabase=$(basename $inlda)
    inldabase=${inldabase%.*}

    inmatbase=$(basename $inmat)
    inmatbase=${inmatbase%.*}

    outbase=${inldabase}.x.${inmatbase}

    outf="${outdir}/${outbase}.RData"

    basesuffix="${jmark}_${binarize}_${Kchoose}_qsub"
    BNAME=${outdir}/${basesuffix}
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "Rscript $rs $inlda $inmat $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $basesuffix

done


