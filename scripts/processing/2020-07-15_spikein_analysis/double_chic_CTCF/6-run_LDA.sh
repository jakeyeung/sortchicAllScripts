#!/bin/sh
# Jake Yeung
# 6-run_LDA.sh
#  
# 2020-08-11


jmem='32G'
jtime='48:00:00'


rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

prefix="K562_spikein"

# outmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins"
# [[ ! -d $outmain0 ]] && mkdir $outmain0
# outmain="${outmain0}/ldaAnalysisBins_${prefix}"
# [[ ! -d $outmain ]] && mkdir $outmain
# [[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA"
# inf="${indir}/count_mat_K36me3_CTCF_dbl_TSS.rds"
# inspike="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA/count_mat_K36me3_CTCF_dbl_CTCFcuts.rds"
outmain="${indir}/LDA"
[[ ! -d $outmain ]] && mkdir $outmain

for inf in `ls -d $indir/count_mat_K36me3_CTCF_dbl_TSS*.rds`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}"
    [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd="cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname --SkipMeanVar"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=LDA_${BNAME} --wrap "$cmd"
done
