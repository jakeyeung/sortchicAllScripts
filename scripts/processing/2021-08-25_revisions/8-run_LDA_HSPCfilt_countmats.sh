#!/bin/sh
# Jake Yeung
# 8-run_LDA_HSPCfilt_countmats.sh
#  
# 2021-09-07

jmem='16G'
jtime='24:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/HSPCs_filt_count_mats"

dnames="TSS_10kb_no_k9 genomewide_bins_50kb"

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

outmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_revisions"
[[ ! -d $outmain0 ]] && mkdir $outmain0


for dname in $dnames; do
    indir=${inmain}/${dname}
    # indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/rds_objs"
    # [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

    for inf in `ls -d $indir/*.rds`; do
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        bname=$(basename $inf)
        bname=${bname%.*}.K-${topicsName}

        prefix="${bname}"
        outmain="${outmain0}/ldaAnalysis_HSPCsOnly_${prefix}"
        [[ ! -d $outmain ]] && mkdir $outmain
        [[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

        outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}"
        [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
        [[ ! -d $outdir ]] && mkdir -p $outdir

        BNAME=$outdir/$bname
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        cmd="cd $workdir; . /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outdir --topics $topics --projname $bname"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=LDA_${bname} --wrap "$cmd"
        # echo $bname
        # echo $cmd
    done

done



