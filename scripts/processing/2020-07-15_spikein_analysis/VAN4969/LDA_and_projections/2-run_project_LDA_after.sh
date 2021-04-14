#!/bin/sh
# Jake Yeung
# 1-run_LDA.sh
#  
# 2020-08-20

# # WRAP UP
# while [[ `squeue -u jyeung | grep LDA | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='32G'
jtime='48:00:00'


# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_LDA_model2.R"
rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cleaned_count_tables_for_lda_and_projections"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
binarize="FALSE"

prefix="K562_spikein_integrate_with_old"

outmain0="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins"
[[ ! -d $outmain0 ]] && mkdir $outmain0

# outmain="${outmain0}/ldaAnalysisBins_${prefix}"
# [[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1

projectmain=${outmain0}/projections_${prefix}
[[ ! -d $projectmain ]] && mkdir $projectmain

ldamain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_K562_spikein_integrate_with_old"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

jbase="count_mats_old_binsize_50000_genomewide"
for jmark in $jmarks; do
    inf="${indir}/${jbase}.${jmark}.old.rds"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    infnew="${indir}/${jbase}.${jmark}.new.rds"
    [[ ! -e $infnew ]] && echo "$infnew not found, exiting" && exit 1

    bname=$(basename $inf)
    bname=${bname%.*}.K-${topicsName}

    ldadir="${ldamain}/lda_outputs.${bname}.binarize.FALSE"
    inlda="${ldadir}/ldaOut.${bname}.Robj"

    projectoutdir=${projectmain}/projections.${bname}
    [[ ! -d $projectoutdir ]] && mkdir $projectoutdir

    outf="${projectoutdir}/${jbase}.projection.RData"

    BNAME=${projectoutdir}/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $infnew $outf"
    # echo $cmd
    echo ${projectoutdir}
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=project_${jmark} --wrap "$cmd" 
done
