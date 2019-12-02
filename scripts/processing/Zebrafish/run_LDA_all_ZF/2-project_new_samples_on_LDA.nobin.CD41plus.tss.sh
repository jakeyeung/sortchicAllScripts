#!/bin/sh
# Jake Yeung
# 2-project_new_samples_on_LDA.nobin.CD41plus.sh
#  
# 2019-11-04

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# jmarks="H3K4me1 H3K4me3 H3K9me3"
jmarks="H3K9me3"

binarize="FALSE"
Kchoose="30"

jmem='32G'
jtime='24:00:00'
# winsizes=100000
winsizes="20000 50000 100000"

for winsize in $winsizes; do
    for jmark in $jmarks; do
        inlda="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-ZFWKM-${jmark}.winsize_${winsize}.merged.K-30.binarize.FALSE/ldaOut.PZ-ChIC-ZFWKM-${jmark}.winsize_${winsize}.merged.K-30.Robj"
        inbase=$(dirname $inlda)
        inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3909_oud3910_merged/countTables_geneTSS/mergedRownamesReadded/PZ-ChIC-ZFWKMCD41plus-${jmark}.winsize_${winsize}.merged.RData"

        [[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
        [[ ! -e $inlda ]] && echo "$inlda not found, continuing" && continue

        outdir=${inbase}/projections
        [[ ! -d $outdir ]] && mkdir $outdir

        inldabase=$(basename $inlda)
        inldabase=${inldabase%.*}

        inmatbase=$(basename $inmat)
        inmatbase=${inmatbase%.*}

        outbase=${inldabase}.x.${inmatbase}

        outf="${outdir}/${outbase}.RData"
        [[ -e $outf ]] && echo "$outf found, continuing" && continue

        basesuffix="${jmark}_${binarize}_${Kchoose}_qsub"
        BNAME=${outdir}/${basesuffix}
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        echo "Rscript $rs $inlda $inmat $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N $basesuffix
        # Rscript $rs $inlda $inmat $outf
        # echo "Rscript $rs $inlda $inmat $outf"
        # echo "Rscript $rs $inlda $inmat $outf"
    done
done


