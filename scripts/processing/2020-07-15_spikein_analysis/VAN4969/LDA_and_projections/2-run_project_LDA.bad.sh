#!/bin/sh
# Jake Yeung
# 2-run_project_LDA.sh
#  
# 2020-08-20

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

binarize="FALSE"
Kchoose="30"

jmem='32G'
jtime='24:00:00'


for jmark in $jmarks; do
    inlda="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-ZFWKM-${jmark}.winsize_${winsize}.merged.K-30.binarize.FALSE/ldaOut.PZ-ChIC-ZFWKM-${jmark}.winsize_${winsize}.merged.K-30.Robj"
    inbase=$(dirname $inlda)
    inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cleaned_count_tables_for_lda_and_projections/count_mats_old_binsize_50000_genomewide.${jmark}.new.rds"

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

    cmd="Rscript $rs $inlda $inmat $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=project_${jmark} --wrap "$cmd" 
done
